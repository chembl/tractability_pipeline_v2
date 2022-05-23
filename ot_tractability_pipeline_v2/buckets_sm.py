#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script contains all functions for the small molecule (SM) workflow.

Created on Mon Feb  3 10:38:51 2020
@author: Melanie Schneider
"""

# import io
# import re
import json
import sys
import time
# import zipfile
# import argparse
# import datetime
import os

# import mygene
import numpy as np
import pandas as pd
import pkg_resources
import itertools
from sqlalchemy import text

PY3 = sys.version > '3'
if PY3:
    import urllib.request as urllib2
else:
    import urllib2


from ot_tractability_pipeline_v2.queries_sm import *

DATA_PATH = pkg_resources.resource_filename('ot_tractability_pipeline_v2', 'data/')


class Small_molecule_buckets(object):
    '''
    Class for assigning genes to tractability buckets
    '''

    ##############################################################################################################
    #
    # Initial setup
    #
    #
    ##############################################################################################################

    def __init__(self, Pipeline_setup, database_url=None, ligand_filter=list()):

        # Cross referencing from Pipeline_setup, prevents repetition for antibody ot_tractability_pipeline
        self.store_fetched = Pipeline_setup.store_fetched

        # list of ensembl IDs for targets to be considered
        self.gene_list = Pipeline_setup.gene_list

        # id_xref dataframe
        self.id_xref = Pipeline_setup.id_xref

        # ChEMBL DB connection
        self.engine = Pipeline_setup.engine

        # # use function from Pipeline setup
        # self.make_request = Pipeline_setup.make_request

        # All chembl data loaded into here
        self.all_chembl_targets = None

        # Unique list of PDB codes:
        self.pdb_list = []

        # Map back to Uniprot accession
        self.pdb_map = {}
        self.acc_map = {}
                
        # URL for PDBe web services
        self.PDB_SERVER_URL = "https://www.ebi.ac.uk/pdbe/api"

        # Load lists for PDB ligand filters (i.e. to remove sugars, solvents etc)
        self.ligand_filter = ligand_filter

        if len(ligand_filter) == 0:
            # Default ligand filter
            # Load unwanted PDB ligands list from text files

            # with open(os.path.join(DATA_PATH, 'organic_solvents.txt')) as solvent_file:
            #     solvents = [a.split('\t')[0] for a in solvent_file]

            # with open(os.path.join(DATA_PATH, 'sugars.txt')) as sugar_file:
            #     sugars = [a.split('\t')[0] for a in sugar_file]

            # with open(os.path.join(DATA_PATH, 'cofactors.txt')) as cofactor_file:
            #     cofactors = [a.split('\t')[0] for a in cofactor_file]

            # with open(os.path.join(DATA_PATH, 'additives.txt')) as additives_file:
            #     additives = [a.split('\t')[0] for a in additives_file]
            
            # # ligand_filter = solvents + sugars + cofactors + additives
            # ligand_filter = solvents + sugars + additives

            # Using unwanted_ligands.txt (mainly buffer compounds, crystallisation agents and sugars) 
            # from "http://cheminfo.u-strasbg.fr:8080/scPDB/2012/db_search/scpdb_generation_process.jsp"
            with open(os.path.join(DATA_PATH, 'unwanted_ligands.txt')) as unwanted_file:
                unwanted = list(itertools.chain.from_iterable([a.split() for a in unwanted_file]))
            
            # unwanted_ligands_added.txt contains additional additives and solvents
            with open(os.path.join(DATA_PATH, 'unwanted_ligands_added.txt')) as unwanted_plus_file:
                unwanted_plus = list(itertools.chain.from_iterable([a.split() for a in unwanted_plus_file]))

            ligand_filter = unwanted + unwanted_plus

            self.ligand_filter = ligand_filter


    ##############################################################################################################
    #
    # Functions relating to buckets 1-3
    # Clinical Compounds
    #
    ##############################################################################################################

    def _search_chembl_clinical(self):
        '''
        Search for all targets in ChEMBL, returning data required for buckets 1-3
        :return:
        '''

        # print("\t- Querying ChEMBL...")
        self.all_chembl_targets = pd.read_sql_query(text(chembl_clinical_targets), self.engine)
        if self.store_fetched: 
            self.all_chembl_targets.to_csv("{}/sm_all_chembl_clinical_targets.csv".format(self.store_fetched))
        small_mol_info = pd.read_sql_query(text(chembl_clinical_small_mol), self.engine)
        if self.store_fetched: 
            small_mol_info.to_csv("{}/sm_all_chembl_clinical_small_mol.csv".format(self.store_fetched))
        self.all_chembl_targets = self.all_chembl_targets.merge(small_mol_info, how='left', on='parent_molregno')

        if self.store_fetched: 
            self.all_chembl_targets.to_csv("{}/sm_all_chembl_targets.csv".format(self.store_fetched))

    def _process_protein_complexes(self):
        '''
        For protein complexes, see if we know the binding subunit, and only keep these
        :return:
        '''

        # print("\t- Processing protein complexes...")
        pc = self.all_chembl_targets[self.all_chembl_targets['target_type'].str.contains("PROTEIN COMPLEX")]
        not_pc = self.all_chembl_targets[~self.all_chembl_targets['target_type'].str.contains("PROTEIN COMPLEX")]

#        n = 1000
#        targets = pc['tid'].unique()
#        chunks = [targets[i:i + n] for i in range(0, len(targets), n)]
#
#        df_list = []
#
#        # Check if binding sites are defined
#        for chunk in chunks:
#            q = '''
#            select distinct bs.site_id, td.tid 
#            from {0}.target_dictionary td, {0}.binding_sites bs
#            where td.tid = bs.tid and td.tid IN {1}
#            '''.format(CHEMBL_VERSION, tuple(chunk))
#            df_list.append(pd.read_sql_query(text(q), self.engine))
#
#        # Merge will set those with unknown binding site as NAN
#        binding_site_info = pd.concat(df_list, sort=False)
#
#        if self.store_fetched: 
#            binding_site_info.to_csv("{}/sm_chembl_binding_site_info.csv".format(self.store_fetched))
#
#        pc = pc.merge(binding_site_info, how='left', on='tid')
        defined = pc[pc['site_id'].notnull()]
        undefined = pc[~pc['site_id'].notnull()]

        # if binding site is defined, only take the subunits that are involved in the binding

        n = 1000
        targets = defined['accession'].unique()
        chunks = [targets[i:i + n] for i in range(0, len(targets), n)]
        df_list2 = []

        for chunk in chunks:
            q2 = '''
            select distinct sc.component_id, cs.accession
            from {0}.component_sequences cs, {0}.site_components sc
            where cs.component_id = sc.component_id
            and cs.accession in {1}'''.format(CHEMBL_VERSION, tuple(chunk))
            df_list2.append(pd.read_sql_query(text(q2), self.engine))

        binding_subunit = pd.concat(df_list2, sort=False)

        if self.store_fetched: 
            binding_subunit.to_csv("{}/sm_chembl_binding_subunit.csv".format(self.store_fetched))

        temp_pc = pc.merge(binding_subunit, on='accession')
        binding_subunits = temp_pc[temp_pc['component_id'].notnull()]

        self.all_chembl_targets = pd.concat([binding_subunits, undefined, not_pc], sort=False)

    def _assess_clinical(self):
        '''
        Merge the results of the ChEMBL search with the OT data (right join, to keep all OT targets)
        Group activity data by target, assign the Max Phase for each targets, and use it to assign buckets 1 to 3

        '''

        print("\t- Assessing clinical buckets 1-3...")

        def set_as_tuple(x):
            return tuple(x)

        def set_strings(x):
            ''' concatenate in string and include only if it is a string (not nan), and exists '''
            return ",".join([y for y in x if isinstance(y,str) and y])

        def set_strings2(x):
            ''' concatenate in string and include only if it is a string (not nan), and exists '''
            return " | ".join([y for y in x if isinstance(y,str) and y])

        #print(self.id_xref['symbol'])
        self._search_chembl_clinical()
        self._process_protein_complexes()

        self.gene_xref = self.id_xref[['accession', 'ensembl_gene_id', 'symbol', 'Entry name', 'Protein names', 
                                       'IDG_family', 'IDG_dtoclass', 'IDG_tdl', 'IDG_fam',
                                       'GO_BioProcess', 'GO_MolFunction']]

        self.clinical_evidence = self.all_chembl_targets.merge(self.gene_xref, how='outer', on='accession')
#        self.out_df = self.all_chembl_targets.merge(self.gene_xref, how='outer', on='accession')
#        self.clinical_evidence = self.out_df

        self.clinical_evidence.to_csv("{}/sm_clinical_evidence.csv".format(self.store_fetched), index=False)

        self.clinical_evidence.drop(['component_id', 'ref_id', 'ref_type', 'tid', 'molregno',
                          'parent_molregno', 'ref_url'], axis=1, inplace=True)

#        # Make sure max phase is for correct indication
#        self.clinical_evidence = self.clinical_evidence[
#            self.clinical_evidence['max_phase'] == self.clinical_evidence['max_phase_for_ind']]

        # pre-processing groupby on two columns to get highest max_phase by drug (and target)
        f0 = {x: 'first' for x in self.clinical_evidence.columns if x not in ['ensembl_gene_id','drug_chembl_id']}
        f0['max_phase_for_ind'] = 'max'
        f0['max_phase'] = 'max'
        
        # as groupby().agg('first') excludes rows with empty column value (nan) - need to fillna and replace afterwards
#        self.clinical_evidence = self.clinical_evidence.groupby([self.clinical_evidence['ensembl_gene_id'],
#                                                                 self.clinical_evidence['drug_chembl_id'].fillna('tmp')
#                                                                 ], as_index=False).agg(f0).replace({'drug_chembl_id':{'tmp': np.nan}})
        self.clinical_evidence = self.clinical_evidence.groupby(['ensembl_gene_id', 'drug_chembl_id'], as_index=False).agg(f0)

        self.out_df = self.gene_xref.merge(self.clinical_evidence.drop(['accession', 'symbol', 'Entry name', 'Protein names', 
                                                                        'IDG_family', 'IDG_dtoclass', 'IDG_tdl', 'IDG_fam',
                                                                        'GO_BioProcess', 'GO_MolFunction'], axis=1), how='outer', on='ensembl_gene_id')

        # copy 'max_phase' column to 'clinical_phase', but first convert to integer (from float), then to string and replace string nan by real nan (that it can correctly be detected during aggregation)
        self.out_df['clinical_phase'] = self.out_df['max_phase'].fillna(-1).astype(int).astype(str).replace('-1',np.nan)


        f = {x: 'first' for x in self.out_df.columns if x != 'ensembl_gene_id'}
        f['max_phase_for_ind'] = 'max'
        f['max_phase'] = 'max'
        f['pref_name'] = set_as_tuple
        f['moa_chembl'] = set_as_tuple
        f['drug_chembl_id'] = set_strings
        f['drug_name'] = set_strings2
        f['clinical_phase'] = set_strings

        self.out_df = self.out_df.groupby(['ensembl_gene_id'], as_index=False).agg(f)#.reset_index(drop=True)

        self.out_df.rename(columns = {'drug_chembl_id':'drug_chembl_ids_sm',
                                      'drug_name':'drug_names_sm',
                                      'clinical_phase':'clinical_phases_sm'}, inplace = True)
        
        self.out_df['Bucket_1_sm'] = 0
        self.out_df['Bucket_2_sm'] = 0
        self.out_df['Bucket_3_sm'] = 0

        self.out_df.loc[(self.out_df['max_phase'] == 4), 'Bucket_1_sm'] = 1
        self.out_df.loc[(self.out_df['max_phase'] < 4) & (self.out_df['max_phase'] >= 2), 'Bucket_2_sm'] = 1
        self.out_df.loc[(self.out_df['max_phase'] < 2) & (self.out_df['max_phase'] > 0), 'Bucket_3_sm'] = 1

        #print(self.out_df['symbol'])
        print(self.out_df.columns)


    ##############################################################################################################
    #
    # Functions relating to bucket 4
    # PDB
    #
    ##############################################################################################################

    # def make_request(self, url, data):
    #     request = urllib2.Request(url)

    #     try:
    #         url_file = urllib2.urlopen(request, data)
    #     except urllib2.HTTPError as e:
    #         if e.code == 404:
    #             print("[NOTFOUND %d] %s" % (e.code, url))
    #         else:
    #             print("[ERROR %d] %s" % (e.code, url))

    #         return None

    #     return url_file.read().decode()

    # Method is used in several workflows
    @staticmethod
    def make_request(url, data):
        request = urllib2.Request(url)

        try:
            url_file = urllib2.urlopen(request, data)
        except urllib2.HTTPError as e:
            if e.code == 404:
                print("[NOTFOUND %d] %s" % (e.code, url))
            else:
                print("[ERROR %d] %s" % (e.code, url))

            return None

        return url_file.read().decode()


    def post_request_pdb(self, url, data, pretty=False):
        full_url = "%s/%s/?pretty=%s" % (self.PDB_SERVER_URL, url, str(pretty).lower())

        if isinstance(data, (list, tuple)):
            data = ",".join(data)

        return self.make_request(full_url, data.encode())

    def _pdb_list(self, s):

        pdb = s['pdb']
        acc = s['accession']

        if not isinstance(pdb, list): pdb = [pdb]

        # Python 2/3 compatability
        try: pdb = [p.lower() for p in pdb if isinstance(p,(str,unicode))] #Python 2
        except: pdb = [p.lower() for p in pdb if isinstance(p,str)] #Python 3

        self.pdb_list += pdb
        for p in pdb:
            try: self.pdb_map[p].add(acc)
            except KeyError: self.pdb_map[p] = {acc}
            
            try: self.acc_map[acc].append(p)
            except KeyError: self.acc_map[acc] = [p]

        self.pdb_list = list(set(self.pdb_list))

    def _has_ligands(self, ligand_li):

        ligands = [l['chem_comp_id'] for l in ligand_li if l['chem_comp_id'] not in self.ligand_filter
                   and 'ION' not in l['chem_comp_name']]

        if len(ligands) > 0:
            return True
        else:
            return False

    def _pdb_ligand_info(self):
        '''
        Use PDBe webservices to get ligand information. POST requests allow up to 1000 PDBs to be submitted each request
        :return:
        '''

        self.id_xref.apply(self._pdb_list, axis=1)

        # Fails with n=1000, runs with n=750
        n = 750
        chunks = [self.pdb_list[i:i + n] for i in range(0, len(self.pdb_list), n)]

        self.no_ligands = []
        self.good_ligands = []
        self.bad_ligands = []

        all_results = {}
        for chunk in chunks:
            ligand_url = '/pdb/entry/ligand_monomers'

            data = ','.join(chunk)

            results = json.loads(self.post_request_pdb(ligand_url, data, False))
            all_results.update(results)
            
            # PDBs without ligands
            chunk_no_ligand = [p for p in chunk if p not in results.keys()]

            # Separate PDBs that do or don't contain suitable ligands
            chunk_good_ligand = [p for p in results.keys() if self._has_ligands(results[p])]
            chunk_bad_ligand = [p for p in results.keys() if not self._has_ligands(results[p])]

            # Add chunk info
            self.no_ligands += chunk_no_ligand
            self.good_ligands += chunk_good_ligand
            self.bad_ligands += chunk_bad_ligand

            time.sleep(1.5)

        if self.store_fetched: 
            json.dump(all_results, open("{}/sm_pdb_ligand_info.json".format(self.store_fetched), 'wt'))


    def _known_pdb_ligand(self, s):

        if s in self.acc_known_lig:
            return list(set([p for p in self.acc_map[s] if p in self.good_ligands]))
        else:
            return np.nan

    def _other_pdb(self, s):

        if s in self.acc_unknown_lig:
            return list(set([p for p in self.acc_map[s] if p in self.no_ligands or p in self.bad_ligands]))
        else:
            return np.nan

    def _assess_pdb(self):
        '''
        Does the target have a ligand-bound protein crystal structure?
        '''

        print("\t- Assessing PDB bucket 4...")

        # Download ligand info from pdb
        self._pdb_ligand_info()

        # Accession numbers with good PDB ligand
        self.acc_known_lig = list({c for pdb in self.good_ligands for c in self.pdb_map[pdb]})
        # Accession numbers with PDB entry (good, bad, or no ligand)
        self.acc_unknown_lig = list({c for pdb in (self.good_ligands + self.no_ligands + self.bad_ligands) for c in self.pdb_map[pdb]})

        self.out_df['PDB_Known_Ligand'] = self.out_df['accession'].apply(self._known_pdb_ligand)
        self.out_df['PDB_other'] = self.out_df['accession'].apply(self._other_pdb)

        self.out_df['Bucket_4_sm'] = 0

        self.out_df.loc[(self.out_df['PDB_Known_Ligand'].notna()), 'Bucket_4_sm'] = 1

        print(self.out_df.columns)


    ##############################################################################################################
    #
    # Functions relating to bucket 5
    # High quality ChEBML compounds 
    # (pchembl_value >=5.5 or if no pchembl_value, then 
    #       any activity <=3000nM or 'k_on'/'Km' measure or 'deltaTm'/'Tm'>=2degreesC or 'Residual activity'<=10% or 'Activity'/'Imax'>=70%, or 'Emax'/'Efficacy'>=120%
    # PFI <=7, structural_alerts <= 2)
    #
    ##############################################################################################################

    def _search_chembl(self):

        self.activities = pd.concat([pd.read_sql_query(text(pchembl_q), self.engine),
                                     pd.read_sql_query(text(nm_q), self.engine),
                                     pd.read_sql_query(text(km_kon_q), self.engine),
                                     pd.read_sql_query(text(D_Tm_q), self.engine),
                                     pd.read_sql_query(text(residual_act_q), self.engine),
                                     pd.read_sql_query(text(Imax_q), self.engine),
                                     pd.read_sql_query(text(Emax_q), self.engine)], sort=False)

        if self.store_fetched: 
            self.activities.to_csv("{}/sm_chembl_activities.csv".format(self.store_fetched))

        # For faster testing
        # self.activities = pd.read_sql_query(text(pchembl_q), self.engine)

    def _structural_alerts(self):
        q = '''
        select distinct cs.canonical_smiles, csa.alert_id /*, sa.alert_name, sas.set_name */
        from {0}.compound_structures cs,
            {0}.compound_structural_alerts csa,
            {0}.structural_alerts sa,
            {0}.molecule_dictionary md,
            {0}.structural_alert_sets sas
        where cs.molregno = md.molregno
        and cs.molregno = csa.molregno
        and csa.alert_id = sa.alert_id
        and sa.alert_set_id = 1
        and sa.alert_set_id = sas.alert_set_id
        '''.format(CHEMBL_VERSION)

        alerts = pd.read_sql_query(text(q), self.engine)

        if self.store_fetched: 
            alerts.to_csv("{}/sm_chembl_alerts.csv".format(self.store_fetched))
        
        # counting the alerts per compound
        alerts = alerts.groupby('canonical_smiles', as_index=False).count()
        alerts.rename(columns={'alert_id': 'alert_id_count'}, inplace=True)

        return alerts

    def _calc_pfi(self, s):

        ar = s['aromatic_rings']
        logd = s['cx_logd']
        return ar + logd

    def _assess_chembl(self):
        '''
        Does the target have ligands in ChEMBL (pchembl_value >= 5.5, PFI <=7, structural_alerts <= 2)
        Scaffold counting currently not implemented (scaffolds >= 2)
        '''

        print("\t- Assessing ChEMBL bucket 5...")
        
        self._search_chembl()
        self.activities['pfi'] = self.activities.apply(self._calc_pfi, axis=1)

        alerts = self._structural_alerts()

        self.activities = self.activities.merge(alerts, how='left', on='canonical_smiles')
        self.activities['alert_id_count'].fillna(0, inplace=True)
        self.activities.to_csv("{}/sm_chembl_activities_processed.csv".format(self.store_fetched))
        df = self.activities[(self.activities['pfi'] <= 7) & (self.activities['alert_id_count'] <= 2)]

        f = {x: 'first' for x in df.columns}
        f['canonical_smiles'] = 'count'
        df2 = df.groupby('accession').agg(f).reset_index(drop=True)
        df2 = df2[['accession', 'canonical_smiles', 'target_chembl_id']]
        self.out_df = df2.merge(self.out_df, how='right', on='accession')
        self.out_df['Bucket_5_sm'] = 0
        self.out_df['target_chembl_id_y'] = self.out_df['target_chembl_id_y'].fillna(self.out_df['target_chembl_id_x'])
        self.out_df.rename(columns={'target_chembl_query_y': 'target_chembl_query'}, inplace=True)

        self.out_df.loc[(self.out_df['canonical_smiles'] >= 2), 'Bucket_5_sm'] = 1

        self.out_df.rename(columns={'canonical_smiles':'High_Quality_ChEMBL_compounds'}, inplace=True)

        print(self.out_df.columns)

        # Use RDKit to count scaffolds
        # PandasTools.AddMoleculeColumnToFrame(self.activities,'canonical_smiles','molecule')
        # PandasTools.AddMurckoToFrame(self.activities,molCol='molecule',MurckoCol='scaffold',Generic=True)
        # self.activities.to_csv('scaffolds.csv')


    ##############################################################################################################
    #
    # Functions relating to buckets 6-7
    # DrugEBIlity
    #
    ##############################################################################################################

    def _assess_pockets(self):
        '''
        Does the target have a DrugEBIlity ensemble score >=0.7 (bucket 6) or  0<score<0.7 (bucket 7)
        '''

        print("\t- Assessing DrugEBIlity buckets 6 and 7...")

        df = pd.read_csv(os.path.join(DATA_PATH, 'drugebility_scores.csv'))

        df = df.merge(self.gene_xref, on='accession', how='right')
        #print(df.dtypes)
        df = df.groupby('ensembl_gene_id', as_index=False).max(numeric_only=True)
        df['ensemble'].fillna(-1, inplace=True)
        df.rename(columns={'ensemble': 'DrugEBIlity_score'}, inplace=True)

        self.out_df = df.merge(self.out_df, how='right', on='ensembl_gene_id', suffixes=['_drop', ''])
        self.out_df['Bucket_6_sm'] = 0
        self.out_df['Bucket_7_sm'] = 0

        self.out_df.loc[(self.out_df['DrugEBIlity_score'] >= 0.7), 'Bucket_6_sm'] = 1
        self.out_df.loc[(self.out_df['DrugEBIlity_score'] > 0) & (self.out_df['DrugEBIlity_score'] < 0.7), 'Bucket_7_sm'] = 1

        print(self.out_df.columns)


    ##############################################################################################################
    #
    # Functions relating to buckets 8
    # Is this target considered druggable using Finan et al's druggable genome?
    #
    ##############################################################################################################

    def _assess_druggable_genome(self):
        '''
        Is this target considered druggable using Finan et al's druggable genome?
        '''

        print("\t- Assessing druggable genome bucket 8...")

        df = pd.read_csv(os.path.join(DATA_PATH, 'druggable_genome.csv'))
        #df = df[['ensembl_gene_id', 'small_mol_druggable']]
        #df.rename(columns={'small_mol_druggable': 'Small_Molecule_Druggable_Genome_Member'}, inplace=True)

        self.out_df = self.out_df.merge(df[['ensembl_gene_id', 'small_mol_druggable']], how='left', on='ensembl_gene_id')
        # some ensembl_gene_id from Finan et al. are not valid any more > need to merge additionally on 'hgnc_names'/'symbol'
        self.out_df = self.out_df.merge(df[['hgnc_names', 'small_mol_druggable']], how='left', left_on='symbol', right_on='hgnc_names')
        
        self.out_df['Small_Molecule_Druggable_Genome_Member'] = self.out_df['small_mol_druggable_x'].fillna(value=self.out_df['small_mol_druggable_y'])
        
        self.out_df['Bucket_8_sm'] = 0
        self.out_df.loc[(self.out_df['Small_Molecule_Druggable_Genome_Member'] == 'Y'), 'Bucket_8_sm'] = 1
        self.out_df['Small_Molecule_Druggable_Genome_Member'].fillna('N', inplace=True)

        print(self.out_df.columns)


    ##############################################################################################################
    #
    # Functions relating to buckets 9
    # Targets with 'chemical' patents in the last 5 years
    #
    ##############################################################################################################

    def _assign_bucket_9(self):
        '''
        Future bucket

        '''
        # print("\t- Assessing  bucket 9...")

        pass


    ##############################################################################################################
    #
    # Higher level functions relating to the overall process
    #
    #
    ##############################################################################################################

    def _summarise_buckets(self):

        '''
        Calculate the best highest bucket for each target
        :return:
        '''

        print("\t- Summarising buckets...")

        self.out_df['Top_bucket_sm'] = 10
        for x in range(8, 0, -1):
            self.out_df.loc[(self.out_df['Bucket_{}_sm'.format(x)] == 1), 'Top_bucket_sm'] = x
            self.out_df['Bucket_{}_sm'.format(x)].fillna(0, inplace=True)

        self.out_df['Bucket_sum_sm'] = self.out_df['Bucket_1_sm'] + self.out_df['Bucket_2_sm'] + self.out_df[
            'Bucket_3_sm'] + self.out_df['Bucket_4_sm'] + self.out_df['Bucket_5_sm'] + self.out_df['Bucket_6_sm'] + self.out_df[
                                        'Bucket_7_sm'] + self.out_df['Bucket_8_sm']

    def _clinical_precedence(self, s):
        return 1 * s['Bucket_1_sm'] + 0.7 * s['Bucket_2_sm'] + 0.2 * s['Bucket_3_sm']

    def _discovery_precedence(self, s):
        return 0.7 * s['Bucket_4_sm'] + 0.3 * s['Bucket_5_sm']

    def _predicted_tractable(self, s):
        return 0.7 * s['Bucket_6_sm'] + 0.3 * s['Bucket_7_sm'] + 0.3 * s['Bucket_8_sm']

    def assign_buckets(self):
        '''
        Assigns the supplied list of gene IDs into their corresponding tractability buckets.
        :return: A Pandas DataFrame containing the Ensembl gene ID and associated tractability bucket
        '''

        self._assess_clinical()
        self._assess_pdb()
        self._assess_chembl()
        self._assess_pockets()
        self._assess_druggable_genome()
        # self._assign_bucket_9()
        self._summarise_buckets()

        self.out_df.index = self.out_df['ensembl_gene_id']
        
        # drop rows without gene symbol
        self.out_df = self.out_df.dropna(subset=['symbol'])
        # drop rows without accession
        self.out_df = self.out_df.dropna(subset=['accession'])
        
        print(self.out_df.columns)
        
        # Add extra buckets to the list below
        self.out_df = self.out_df[['symbol', 'accession', 'Entry name', 'Protein names', 
                                   'IDG_family', 'IDG_dtoclass', 'IDG_tdl', 'IDG_fam',
                                   'GO_BioProcess', 'GO_MolFunction',
                                   'Bucket_1_sm', 'Bucket_2_sm', 'Bucket_3_sm', 
                                   'Bucket_4_sm', 'Bucket_5_sm', 'Bucket_6_sm', 
                                   'Bucket_7_sm', 'Bucket_8_sm', 
                                   'Bucket_sum_sm', 'Top_bucket_sm',
                                   'drug_chembl_ids_sm', 'drug_names_sm', 'clinical_phases_sm',
                                   'PDB_Known_Ligand', 'PDB_other',
                                   'High_Quality_ChEMBL_compounds', 'DrugEBIlity_score', 'Small_Molecule_Druggable_Genome_Member' 
                                   ]]

        # Calculate category scores and assign highest category to each target

        self.out_df['Category_sm'] = 'Unknown'
        self.out_df['Clinical_Precedence_sm'] = self.out_df.apply(self._clinical_precedence, axis=1)
        self.out_df['Discovery_Precedence_sm'] = self.out_df.apply(self._discovery_precedence, axis=1)
        self.out_df['Predicted_Tractable_sm'] = self.out_df.apply(self._predicted_tractable, axis=1)

        self.out_df.loc[(self.out_df['Top_bucket_sm'] <= 3), 'Category_sm'] = 'Clinical_Precedence_sm'
        self.out_df.loc[(self.out_df['Top_bucket_sm'] == 4) | (self.out_df['Top_bucket_sm'] == 5),
                        'Category_sm'] = 'Discovery_Precedence_sm'
        self.out_df.loc[
            (self.out_df['Top_bucket_sm'] == 6) | (self.out_df['Top_bucket_sm'] == 7) | (self.out_df['Top_bucket_sm'] == 8),
            'Category_sm'] = 'Predicted_Tractable_sm'

        # self.out_df.sort_values(['Clinical_Precedence_sm', 'Discovery_Precedence_sm', 'Predicted_Tractable_sm'],
        #                         ascending=[False, False, False], inplace=True)

        # Cleaning columns
        self.out_df['drug_chembl_ids_sm'].fillna('', inplace=True)
        self.out_df['drug_names_sm'].fillna('', inplace=True)
        self.out_df['clinical_phases_sm'].fillna('', inplace=True)

        # create dictionaries from 'drug_chembl_ids_sm' and 'clinical_phases_sm'/'drug_names_sm'
        self.out_df['drug_names_dict_sm'] = self.out_df.apply(lambda row : dict(zip(row['drug_chembl_ids_sm'].split(","), row['drug_names_sm'].split(" | "))), axis=1)
        self.out_df['clinical_phases_dict_sm'] = self.out_df.apply(lambda row : dict(zip(row['drug_chembl_ids_sm'].split(","), row['clinical_phases_sm'].split(","))), axis=1)
        
        # Cleaning column: setting selected culumns in list format to improve visualization e.g. with Excel
        # and remove duplicates while keeping order using "list(dict.fromkeys(lst))"
        self.out_df['drug_chembl_ids_sm'] = self.out_df['drug_chembl_ids_sm'].apply(lambda x: list(dict.fromkeys(x.split(","))))
        self.out_df['drug_names_sm'] = self.out_df['drug_names_sm'].apply(lambda x: list(dict.fromkeys(x.split(" | "))))
       
        print(self.out_df.columns)

        # return self.out_df
        return self.out_df.astype({x: 'int64' for x in self.out_df.columns if "Bucket" in x})



    @staticmethod
    def sm2json(d):
        return {
            'Bucket_scores': {'Bucket_1_sm':d.Bucket_1_sm, 'Bucket_2_sm':d.Bucket_2_sm, 'Bucket_3_sm':d.Bucket_3_sm, 
                              'Bucket_4_sm':d.Bucket_4_sm, 'Bucket_5_sm':d.Bucket_5_sm, 'Bucket_6_sm':d.Bucket_6_sm, 
                              'Bucket_7_sm':d.Bucket_7_sm, 'Bucket_8_sm':d.Bucket_8_sm}, #, 'Bucket_9_sm':d.Bucket_9_sm
            'Bucket_evaluation': {'Bucket_sum_sm':d.Bucket_sum_sm, 'Top_bucket_sm':d.Top_bucket_sm},
            'Category_scores': {'Clinical_Precedence_sm':d.Clinical_Precedence_sm, 
                                'Discovery_Precedence_sm':d.Discovery_Precedence_sm, 
                                'Predicted_Tractable_sm':d.Predicted_Tractable_sm}, 
            'Category_evaluation': {'Top_Category_sm':d.Category_sm},
#            'Bucket_assessment_data': {'Bucket_1_2_3_sm': {'drug_chembl_ids_sm':d.drug_chembl_ids_sm, 'drug_names_sm':d.drug_names_sm}, 
            'Bucket_assessment_data': {'Bucket_1_2_3_sm': {'drugs_in_clinic':d.drug_names_dict_sm, 'max_clinical_phase':d.clinical_phases_sm}, 
                                     'Bucket_4_sm': {'PDB_known_ligands':d.PDB_Known_Ligand}, 
                                     'Bucket_5_sm': {'number_high_quality_ChEMBL_compounds':d.High_Quality_ChEMBL_compounds}, 
                                     'Bucket_6-7_sm': {'DrugEBIlity_score':d.DrugEBIlity_score}, 
                                     'Bucket_8_sm': {'Small_Molecule_Druggable_Genome_Member':d.Small_Molecule_Druggable_Genome_Member}}
            }

