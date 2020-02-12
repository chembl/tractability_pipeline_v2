#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  3 10:38:51 2020

@author: Melanie Schneider
"""
# import json
# import time
# import argparse
# import datetime
import io
import re
import sys
import zipfile
import os

# import mygene
# import numpy as np
import pandas as pd
from pandas.io.json import json_normalize
import pkg_resources
# from sqlalchemy import create_engine

PY3 = sys.version > '3'
if PY3:
    import urllib.request as urllib2
else:
    import urllib2


from ot_tractability_pipeline_v2.queries_ab import *

DATA_PATH = pkg_resources.resource_filename('ot_tractability_pipeline_v2', 'data/')

class Antibody_buckets(object):
    '''
    Class for assigning genes to tractability buckets
    '''

    ##############################################################################################################
    #
    # Initial setup
    #
    #
    ##############################################################################################################

    def __init__(self, Pipeline_setup, database_url=None, prev_output=None):

        # Cross referencing from Pipeline_setup, prevents repetition for antibody ot_tractability_pipeline
        self.store_fetched = Pipeline_setup.store_fetched

        # list of ensembl IDs for targets to be considered
        self.gene_list = Pipeline_setup.gene_list

        # id_xref dataframe
        self.id_xref = Pipeline_setup.id_xref

        # ChEMBL DB connection
        self.engine = Pipeline_setup.engine

        # All chembl data loaded into here
        self.all_chembl_targets = None

        # Unique list of PDB codes:
        self.pdb_list = []

        # Map back to Uniprot accession
        self.pdb_map = {}

        # Load accepted GO locations
        self.accepted_go_locs = {}
        with open(os.path.join(DATA_PATH, 'go_accepted_loc.tsv')) as go_loc_file:
            for line in go_loc_file:
                line = line.split('\t')
                self.accepted_go_locs[line[0]] = line[2]

        # If antibody results are to be combined with small molecule results, append antibody columns to sm results
        # Otherwise, use the id_xref dataframe
        if prev_output is not None:
            go_data = self.id_xref[['ensembl_gene_id', 'go.CC']]
            self.out_df = prev_output.merge(go_data, how='outer', on='ensembl_gene_id')

        else:
            self.out_df = self.id_xref

    ##############################################################################################################
    #
    # Functions relating to buckets 1-3
    # Clinical antibodies
    #
    ##############################################################################################################

    def _process_protein_complexes(self):
        '''
        For protein complexes, see if we know the binding subunit, and only keep these

        :return:
        '''

        pc = self.all_chembl_targets[self.all_chembl_targets['target_type'].str.contains("PROTEIN COMPLEX")]
        not_pc = self.all_chembl_targets[~self.all_chembl_targets['target_type'].str.contains("PROTEIN COMPLEX")]

        defined = pc[pc['site_id'].notnull()]
        undefined = pc[~pc['site_id'].notnull()]

        # if binding site is defined, only take the subunits that are involved in the binding

        n = 1000
        targets = defined['site_id'].unique()
        chunks = [targets[i:i + n] for i in range(0, len(targets), n)]
        df_list2 = []

        for chunk in chunks:
            q2 = '''
            select distinct cs.accession, sc.component_id
            from {0}.component_sequences cs,
                {0}.site_components sc
            where cs.component_id = sc.component_id
            and sc.site_id in {1}'''.format(CHEMBL_VERSION, tuple(chunk))
            df_list2.append(pd.read_sql_query(q2, self.engine))

        binding_subunit = pd.concat(df_list2, sort=False)

        if self.store_fetched: 
            binding_subunit.to_csv("{}/ab_chembl_binding_subunit.csv".format(self.store_fetched))

        temp_pc = pc.merge(binding_subunit, on='accession')
        binding_subunits = temp_pc[temp_pc['component_id'].notnull()]

        self.all_chembl_targets = pd.concat([binding_subunits, undefined, not_pc], sort=False)

    def _assign_buckets_1_to_3(self):
        '''
        Merge the results of the ChEMBL search with the OT data (right join, to keep all OT targets)

        Group activity data by target, assign the Max Phase for each targets, and use it to assign buckets 1 to 3

        Possible duplication of effort of the OT known_drug score

        :return:
        '''

        self.all_chembl_targets = pd.read_sql_query(chembl_clinical_ab_targets, self.engine)
        self.all_chembl_targets.loc[self.all_chembl_targets['ref_type'] == 'Expert', ['ref_id', 'ref_url']] = 'NA'

        #

        self._process_protein_complexes()

        ab_info = pd.read_sql_query(chembl_clinical_ab, self.engine)
        self.all_chembl_targets = self.all_chembl_targets.merge(ab_info, how='left', on='parent_molregno')

        if self.store_fetched: 
            self.all_chembl_targets.to_csv("{}/ab_all_chembl_targets.csv".format(self.store_fetched))

        # Make sure max phase is for correct indication

        self.all_chembl_targets = self.all_chembl_targets[
            self.all_chembl_targets['max_phase'] == self.all_chembl_targets['max_phase_for_ind']]

        def other_func(x):
            return tuple(x)

        f = {x: other_func for x in self.all_chembl_targets if x != 'accession'}
        f['max_phase_for_ind'] = 'max'
        f['max_phase'] = 'max'

        self.all_chembl_targets = self.all_chembl_targets.groupby('accession', as_index=False).agg(f)

        self.out_df = self.all_chembl_targets.merge(self.out_df, how='outer', on='accession', suffixes=('_sm', '_ab'))

        self.out_df.drop(['component_id', 'drug_name', 'ref_id', 'ref_type', 'tid',
                          'ref_url'], axis=1, inplace=True)

        f = {x: 'first' for x in self.out_df.columns if x != 'ensembl_gene_id'}
        f['max_phase_for_ind'] = 'max'

        self.out_df = self.out_df.groupby(['ensembl_gene_id'], as_index=False).agg(f)

        self.out_df['Bucket_1_ab'] = 0
        self.out_df['Bucket_2_ab'] = 0
        self.out_df['Bucket_3_ab'] = 0

        self.out_df.loc[(self.out_df['max_phase_for_ind'] == 4), 'Bucket_1_ab'] = 1
        self.out_df.loc[
            (self.out_df['max_phase_for_ind'] < 4) & (self.out_df['max_phase_for_ind'] >= 2), 'Bucket_2_ab'] = 1
        self.out_df.loc[
            (self.out_df['max_phase_for_ind'] < 2) & (self.out_df['max_phase_for_ind'] > 0), 'Bucket_3_ab'] = 1

    ##############################################################################################################
    #
    # Functions relating to buckets 4, 6 and 7
    # Uniprot location
    #
    ##############################################################################################################

    def _make_request(self, url, data):
        request = urllib2.Request(url)

        try:
            url_file = urllib2.urlopen(request)
        except urllib2.HTTPError as e:
            if e.code == 404:
                print("[NOTFOUND %d] %s" % (e.code, url))
            else:
                print("[ERROR %d] %s" % (e.code, url))

            return None

        return url_file.read().decode()

    def _post_request(self, url, data):
        base = 'http://www.uniprot.org'
        full_url = "%s/%s" % (base, url)

        if isinstance(data, (list, tuple)):
            data = ",".join(data)

        return self._make_request(full_url, data.encode())

    def split_loc(self, s):
        '''
        Process the delimited string returned from uniprot webservice call
        :param s:
        :return:
        '''
        try:
            loc_list = [a.split('. ') for a in s.split(';')]
            loc_evidence_li = []

        except AttributeError:
            return [('na', 'na')]

        for l in loc_list:
            for l2 in l:
                if l2 == '':
                    continue

                l2 = l2.strip()
                if l2.startswith('Note'):
                    break

                l2 = l2.replace('SUBCELLULAR LOCATION: ', '')
                # evidence = l2[l2.find("{") + 1:l2.find("}")]
                evidence = re.findall(r'\{([^]]*)\}', l2)
                evidence = [e for x in evidence for e in x.split(',')]

                if len(evidence) == 0:
                    evidence = ['Unknown evidence type']

                locations = l2.split('{')[0]

                if locations != '':
                    loc_evidence_li.append((evidence, locations))
        return loc_evidence_li

    def _check_evidence(self, evidence_li):

        high_conf_evidence = [e for e in evidence_li if ('ECO:0000269' in e or 'ECO:0000305' in e)]

        if len(high_conf_evidence) > 0:
            return True
        return False

    def _set_b4_flag(self, s):

        accepted_uniprot_high_conf = [a[1] for a in s['Subcellular location [CC]'] if
                                      ('Cell membrane' in a[1] or 'Secreted' in a[1]) and (self._check_evidence(a[0]))]

        all_uniprot_high_conf = [(a[1], a[0]) for a in s['Subcellular location [CC]'] if self._check_evidence(a[0])]

        if len(accepted_uniprot_high_conf) > 0:
            b4_flag = 1
        else:
            b4_flag = 0

        return b4_flag, all_uniprot_high_conf

    def _set_b6_flag(self, s):

        accepted_uniprot_med_conf = [a[1] for a in s['Subcellular location [CC]'] if
                                     ('Cell membrane' in a[1] or 'Secreted' in a[1]) and not self._check_evidence(a[0])]

        all_uniprot_med_conf = [(a[1], a[0]) for a in s['Subcellular location [CC]'] if not self._check_evidence(a[0])]

        if len(accepted_uniprot_med_conf) > 0:
            b6_flag = 1
        else:
            b6_flag = 0

        return b6_flag, all_uniprot_med_conf

    def _assign_bucket_4_and_6(self):
        '''
        Uniprot (loc): Targets in "Cell membrane" or "Secreted", high confidence
        '''

        # Return all reviewed and Human targets
        url = "uniprot/?format=tab&query=*&fil=reviewed%3ayes+AND+organism%3a%22Homo+sapiens+(Human)+%5b9606%5d%22&columns=id,comment(SUBCELLULAR+LOCATION),comment(DOMAIN),feature(DOMAIN+EXTENT),feature(INTRAMEMBRANE),feature(TOPOLOGICAL+DOMAIN),feature(TRANSMEMBRANE),feature(SIGNAL)"
        data = ['P42336', 'P60484']

        location = self._post_request(url, data)
        location = [x.split('\t') for x in location.split('\n')]
        df = pd.DataFrame(location[1:], columns=location[0])
        df['uniprot_loc_test'] = df['Subcellular location [CC]']

        df['Subcellular location [CC]'] = df['Subcellular location [CC]'].apply(self.split_loc)

        if self.store_fetched: 
            df.to_csv("{}/ab_uniprot_for_buckets_4_and_6.csv".format(self.store_fetched))

        df['Bucket_4_ab'], df['Uniprot_high_conf_loc'] = zip(*df.apply(self._set_b4_flag, axis=1))
        df['Bucket_6_ab'], df['Uniprot_med_conf_loc'] = zip(*df.apply(self._set_b6_flag, axis=1))
        df.rename(columns={'Entry': 'accession'}, inplace=True)

        self.out_df = self.out_df.merge(df, how='left', on='accession')

    ##############################################################################################################
    #
    # Functions relating to buckets 5 and 8
    # GO Cell Component
    #
    ##############################################################################################################

    def _set_b5_b8_flag(self, s):
        try:
            cc = s['go.CC']
        except:
            return 0, [], 0, []

        # Confidence for each evidence type
        evidence_types = {'EXP': 'High', 'IDA': 'High', 'IPI': 'High', 'TAS': 'High', 'IMP': 'High', 'IGI': 'High',
                          'IEP': 'High',
                          'ISS': 'Medium', 'ISO': 'Medium', 'ISA': 'Medium', 'ISM': 'Medium', 'IGC': 'Medium',
                          'IBA': 'Medium', 'IBD': 'Medium', 'IKR': 'Medium', 'IRD': 'Medium', 'RCA': 'Medium',
                          'IEA': 'Medium',
                          'NAS': 'Low', 'IC': 'Low', 'ND': 'Low', 'NR': 'Low'
                          }

        high_conf_loc = []
        med_conf_loc = []
        accepted_high_conf_loc = []
        accepted_med_conf_loc = []

        if isinstance(cc, dict):
            cc = [cc]

        if not isinstance(cc, list):
            return 0, [], 0, []

        for c_dict in cc:
            try:
                go_id = c_dict['id']
                go_loc = c_dict['term']
                evidence = c_dict['evidence']
            except TypeError:
                continue
            try:
                confidence = evidence_types[evidence]
            except KeyError:
                confidence = None

            if confidence == 'High':
                high_conf_loc.append((go_loc, evidence))
            elif confidence == 'Medium':
                med_conf_loc.append((go_loc, evidence))

            if go_id in self.accepted_go_locs.keys():
                if confidence == 'High':
                    accepted_high_conf_loc.append(go_loc)
                elif confidence == 'Medium':
                    accepted_med_conf_loc.append(go_loc)

        b5_flag = 0
        b8_flag = 0

        if len(accepted_high_conf_loc) > 0:
            b5_flag = 1
        elif len(accepted_med_conf_loc) > 0:
            b8_flag = 1

        return b5_flag, high_conf_loc, b8_flag, med_conf_loc

    def _assign_bucket_5_and_8(self):
        '''
        GO CC
        '''
        self.out_df['Bucket_5_ab'], self.out_df['GO_high_conf_loc'], self.out_df['Bucket_8_ab'], self.out_df[
            'GO_med_conf_loc'] = zip(*self.out_df.apply(self._set_b5_b8_flag, axis=1))

    ##############################################################################################################
    #
    # Functions relating to bucket 7
    # Uniprot transmembrane and signal peptides
    #
    ##############################################################################################################

    def _split_loc_b7(self, s):

        try:
            loc_list = [a.split('.') for a in s.split(';')]
            loc_evidence_li = []

        except AttributeError:
            return [('na', 'na')]

        for l in loc_list:
            for l2 in l:
                if l2 == '':
                    continue

                l2 = l2.strip()
                l2 = re.sub(r'{[^}]+}', '', l2)
                if l2.startswith('TRANSMEM') or l2.startswith('SIGNAL'):
                    loc_evidence_li.append(l2)
        return loc_evidence_li

    def _assign_bucket_7(self):
        '''
        Uniprot (SigP + TMHMM): targets with predicted Signal Peptide or Trans-membrane regions, and not destined to
        organelles

        '''
        self.out_df['Bucket_7_ab'] = 0

        self.out_df.loc[(self.out_df['Transmembrane'].str.contains('TRANSMEM', na=False)), 'Bucket_7_ab'] = 1
        self.out_df.loc[(self.out_df['Signal peptide'].str.contains('SIGNAL', na=False)), 'Bucket_7_ab'] = 1

        self.out_df['Transmembrane'] = self.out_df['Transmembrane'].apply(self._split_loc_b7)
        self.out_df['Signal peptide'] = self.out_df['Signal peptide'].apply(self._split_loc_b7)

    ##############################################################################################################
    #
    # Functions relating to buckets 9
    # Human protein atlas - Main location
    #
    ##############################################################################################################

    def _main_location(self, s):

        if s['Reliability'] == 'Validated':
            return s['Validated']
        else:
            return s['Supported']

    def _assign_bucket_9(self):
        '''
        HPA
        '''

        # Download latest file
        zip_file = urllib2.urlopen('https://www.proteinatlas.org/download/subcellular_location.tsv.zip')
        with zipfile.ZipFile(io.BytesIO(zip_file.read()), 'r') as pa_file:
            with pa_file.open('subcellular_location.tsv') as subcell_loc:
                df = pd.read_csv(subcell_loc, sep='\t', header=0)

        df.rename(columns={'Gene': 'ensembl_gene_id'}, inplace=True)

        if self.store_fetched: 
            df.to_csv("{}/ab_proteinatlas_for_bucket_9.csv".format(self.store_fetched))

        df['main_location'] = df.apply(self._main_location, axis=1)
        reliable = df[(df['Reliability'] == 'Supported') | (df['Reliability'] == 'Validated')]

        reliable.drop(columns=["Approved", "Supported", "Uncertain",
                               "Cell cycle dependency", "GO id"], axis=1, inplace=True)

        self.out_df = self.out_df.merge(reliable, on='ensembl_gene_id', how='left')

        self.out_df['Bucket_9_ab'] = 0
        self.out_df.loc[(self.out_df['main_location'].str.contains("Plasma membrane", na=False)), 'Bucket_9_ab'] = 1

    ##############################################################################################################
    #
    # Higher level functions relating to the overall process
    #
    #
    ##############################################################################################################

    def _clinical_precedence(self, s):
        return 1 * s['Bucket_1_ab'] + 0.7 * s['Bucket_2_ab'] + 0.2 * s['Bucket_3_ab']

    def _high_conf_pred(self, s):
        return 0.7 * s['Bucket_4_ab'] + 0.3 * s['Bucket_5_ab']

    def _med_conf_pred(self, s):
        return 0.4 * s['Bucket_6_ab'] + 0.25 * s['Bucket_7_ab'] + 0.25 * s['Bucket_8_ab'] + 0.1 * s['Bucket_9_ab']

    def _summarise_buckets(self):

        self.out_df.drop('go.CC', inplace=True, axis=1)

        self.out_df['Top_bucket_ab'] = 10
        for x in range(9, 0, -1):
            self.out_df.loc[(self.out_df['Bucket_{}_ab'.format(x)] == 1), 'Top_bucket_ab'] = x
            self.out_df['Bucket_{}_ab'.format(x)].fillna(0, inplace=True)

        self.out_df['Bucket_sum_ab'] = self.out_df['Bucket_1_ab'] + self.out_df['Bucket_2_ab'] + self.out_df[
            'Bucket_3_ab'] + self.out_df['Bucket_4_ab'] + self.out_df['Bucket_5_ab'] + self.out_df['Bucket_6_ab'
                                       ] + self.out_df['Bucket_7_ab'] + self.out_df['Bucket_8_ab'] + self.out_df[
                                           'Bucket_9_ab']

        self.out_df.set_index('ensembl_gene_id')

    def assign_buckets(self):
        '''
        Assigns the supplied list of gene IDs into their corresponding tractability buckets.
        :return: A Pandas DataFrame containing the Ensembl gene ID and associated tractability bucket
        '''

        self._assign_buckets_1_to_3()
        self._assign_bucket_4_and_6()
        self._assign_bucket_5_and_8()
        self._assign_bucket_7()
        self._assign_bucket_9()

        self._summarise_buckets()

        # try:

        self.out_df.index = self.out_df['ensembl_gene_id']

        # Columns to keep. This includes columns from the small molecule pipeline

        self.out_df = self.out_df[['symbol', 'accession',
                                   'Bucket_1', 'Bucket_2', 'Bucket_3', 'Bucket_4',
                                   'Bucket_5', 'Bucket_6', 'Bucket_7',
                                   'Bucket_8', 'Bucket_sum', 'Top_bucket', 'Category',
                                   'Clinical_Precedence', 'Discovery_Precedence', 'Predicted_Tractable',
                                   'PDB_Known_Ligand',
                                   'ensemble', 'canonical_smiles', 'small_mol_druggable',
                                   'Bucket_1_ab', 'Bucket_2_ab', 'Bucket_3_ab', 'Bucket_4_ab',
                                   'Bucket_5_ab', 'Bucket_6_ab', 'Bucket_7_ab',
                                   'Bucket_8_ab', 'Bucket_9_ab', 'Bucket_sum_ab', 'Top_bucket_ab',
                                   'Uniprot_high_conf_loc', 'GO_high_conf_loc',
                                   'Uniprot_med_conf_loc',
                                   'GO_med_conf_loc', 'Transmembrane', 'Signal peptide', 'main_location'
                                   ]]
        self.out_df.rename(columns={'canonical_smiles': 'High_Quality_ChEMBL_compounds',
                                    'small_mol_druggable': 'Small_Molecule_Druggable_Genome_Member',
                                    'main_location': 'HPA_main_location', 'Signal peptide': 'Signal_peptide'},
                           inplace=True)
        self.out_df.sort_values(['Clinical_Precedence', 'Discovery_Precedence', 'Predicted_Tractable'],
                                ascending=[False, False, False], inplace=True)

        # Score each category, and label highest category
        self.out_df['Clinical_Precedence_ab'] = self.out_df.apply(self._clinical_precedence, axis=1)
        self.out_df['Predicted_Tractable__High_confidence'] = self.out_df.apply(self._high_conf_pred, axis=1)
        self.out_df['Predicted_Tractable__Medium_to_low_confidence'] = self.out_df.apply(self._med_conf_pred, axis=1)

        self.out_df['Category_ab'] = 'Unknown'

        self.out_df.loc[(self.out_df['Top_bucket_ab'] <= 3), 'Category_ab'] = 'Clinical_Precedence'
        self.out_df.loc[(self.out_df['Top_bucket_ab'] == 4) | (self.out_df['Top_bucket_ab'] == 5),
                        'Category_ab'] = 'Predicted_Tractable__High_confidence'

        self.out_df.loc[
            (self.out_df['Top_bucket_ab'] == 6) | (self.out_df['Top_bucket_ab'] == 7) | (
                    self.out_df['Top_bucket_ab'] == 8) | (self.out_df['Top_bucket_ab'] == 9),
            'Category_ab'] = 'Predicted_Tractable__Medium_to_low_confidence'

        # self.out_df = self.out_df[(self.out_df['Top_bucket'] < 9 ) | (self.out_df['Top_bucket_ab'] < 10) ]

        return self.out_df.astype({x: 'int64' for x in self.out_df.columns if "Bucket" in x})

