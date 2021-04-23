#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script contains all functions for the OtherClinical workflow.

Created on Mon Feb  3 10:38:51 2020
@author: Melanie Schneider
"""

# import io
# import json
# import re
# import time
# import zipfile
# import argparse
# import datetime
# import os
import sys

# import mygene
import numpy as np
import pandas as pd
import pkg_resources
# from sqlalchemy import create_engine

PY3 = sys.version > '3'
if PY3:
    import urllib.request as urllib2
else:
    import urllib2


from ot_tractability_pipeline_v2.queries_othercl import *

DATA_PATH = pkg_resources.resource_filename('ot_tractability_pipeline_v2', 'data/')

class Othercl_buckets(object):
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

        # print("\t\t- Processing protein complexes...")
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
            binding_subunit.to_csv("{}/othercl_chembl_binding_subunit.csv".format(self.store_fetched))

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

        print("\t- Assessing clinical buckets 1-3...")

        self.all_chembl_targets = pd.read_sql_query(chembl_clinical_other_targets, self.engine)
        self.all_chembl_targets.loc[self.all_chembl_targets['ref_type'] == 'Expert', ['ref_id', 'ref_url']] = 'NA'

        self._process_protein_complexes()

        othercl_info = pd.read_sql_query(chembl_clinical_other, self.engine)
        self.all_chembl_targets = self.all_chembl_targets.merge(othercl_info, how='left', on='parent_molregno')

        if self.store_fetched: 
            self.all_chembl_targets.to_csv("{}/othercl_all_chembl_targets.csv".format(self.store_fetched))

#        # Make sure max phase is for correct indication
#        self.all_chembl_targets = self.all_chembl_targets[
#            self.all_chembl_targets['max_phase'] == self.all_chembl_targets['max_phase_for_ind']]

        # pre-processing groupby on two columns to get highest max_phase by drug (and target)
        f0 = {x: 'first' for x in self.all_chembl_targets.columns if x not in ['accession','drug_chembl_id']}
        f0['max_phase_for_ind'] = 'max'
        f0['max_phase'] = 'max'

        self.all_chembl_targets = self.all_chembl_targets.groupby(['accession','drug_chembl_id'], as_index=False).agg(f0)

        def set_as_tuple(x):
            return tuple(x)

        def set_strings(x):
            ''' concatenate in string and include only if it is a string (not nan), and exists '''
            return ",".join([y for y in x if isinstance(y,str) and y])

        def set_strings2(x):
            ''' concatenate in string and include only if it is a string (not nan), and exists '''
            return " | ".join([y for y in x if isinstance(y,str) and y])

        # copy 'max_phase' column to 'clinical_phase', but first convert to integer (from float), then to string and replace string nan by real nan (that it can correctly be detected during aggregation)
        self.all_chembl_targets['clinical_phase'] = self.all_chembl_targets['max_phase'].fillna(-1).astype(int).astype(str).replace('-1',np.nan)

        f = {x: set_as_tuple for x in self.all_chembl_targets if x != 'accession'}
        f['max_phase_for_ind'] = 'max'
        f['max_phase'] = 'max'
        f['drug_chembl_id'] = set_strings
        f['drug_name'] = set_strings2
        f['clinical_phase'] = set_strings

        self.all_chembl_targets = self.all_chembl_targets.groupby('accession', as_index=False).agg(f)

        self.out_df = self.all_chembl_targets.merge(self.out_df, how='outer', on='accession', suffixes=('', '_othercl'))

        self.out_df.drop(['component_id', 'ref_id', 'ref_type', 'tid',
                          'ref_url'], axis=1, inplace=True)

        # self.out_df.to_csv("{}/othercl_out_df_checkpoint1.csv".format(self.store_fetched), index=False)


        f = {x: 'first' for x in self.out_df.columns if x != 'ensembl_gene_id'}
        f['max_phase_for_ind'] = 'max'
        f['max_phase'] = 'max'
        f['drug_chembl_id'] = set_strings
        f['drug_name'] = set_strings2
        f['clinical_phase'] = set_strings

        self.out_df = self.out_df.groupby(['ensembl_gene_id'], as_index=False).agg(f)

        # self.out_df.to_csv("{}/othercl_out_df_checkpoint2.csv".format(self.store_fetched), index=False)

        self.out_df.rename(columns = {'drug_chembl_id':'drug_chembl_ids_othercl',
                                      'drug_name':'drug_names_othercl',
                                      'clinical_phase':'clinical_phases_othercl'}, inplace = True)
        
        self.out_df['Bucket_1_othercl'] = 0
        self.out_df['Bucket_2_othercl'] = 0
        self.out_df['Bucket_3_othercl'] = 0

        self.out_df.loc[(self.out_df['max_phase'] == 4), 'Bucket_1_othercl'] = 1
        self.out_df.loc[
            (self.out_df['max_phase'] < 4) & (self.out_df['max_phase'] >= 2), 'Bucket_2_othercl'] = 1
        self.out_df.loc[
            (self.out_df['max_phase'] < 2) & (self.out_df['max_phase'] > 0), 'Bucket_3_othercl'] = 1


    ##############################################################################################################
    #
    # Higher level functions relating to the overall process
    #
    #
    ##############################################################################################################

    def _clinical_precedence(self, s):
        return 1 * s['Bucket_1_othercl'] + 0.7 * s['Bucket_2_othercl'] + 0.2 * s['Bucket_3_othercl']


    def _summarise_buckets(self):

        print("\t- Summarising buckets...")

        self.out_df.drop('go.CC', inplace=True, axis=1)

        self.out_df['Top_bucket_othercl'] = 10
        for x in range(3, 0, -1):
            self.out_df.loc[(self.out_df['Bucket_{}_othercl'.format(x)] == 1), 'Top_bucket_othercl'] = x
            self.out_df['Bucket_{}_othercl'.format(x)].fillna(0, inplace=True)

        self.out_df['Bucket_sum_othercl'] = self.out_df['Bucket_1_othercl'] + self.out_df['Bucket_2_othercl'] + self.out_df['Bucket_3_othercl']

        self.out_df.set_index('ensembl_gene_id')

    def assign_buckets(self):
        '''
        Assigns the supplied list of gene IDs into their corresponding tractability buckets.
        :return: A Pandas DataFrame containing the Ensembl gene ID and associated tractability bucket
        '''

        self._assign_buckets_1_to_3()
        print(self.out_df.columns)

        self._summarise_buckets()
        print(self.out_df.columns)

        # try:

        self.out_df.index = self.out_df['ensembl_gene_id']

        # drop rows without gene symbol
        self.out_df = self.out_df.dropna(subset=['symbol'])
        # drop rows without accession
        self.out_df = self.out_df.dropna(subset=['accession'])

        # Columns to keep
        self.out_df = self.out_df[['accession', 'symbol',
                                   # 'Bucket_1_sm', 'Bucket_2_sm', 'Bucket_3_sm', 
                                   # 'Bucket_4_sm', 'Bucket_5_sm', 'Bucket_6_sm', 'Bucket_7_sm',
                                   # 'Bucket_8_sm', 'Bucket_sum_sm', 'Top_bucket_sm', 'Category_sm',
                                   # 'Clinical_Precedence_sm', 'Discovery_Precedence_sm', 'Predicted_Tractable_sm',
                                   # 'PDB_Known_Ligand',
                                   # 'ensemble', 'canonical_smiles', 'small_mol_druggable',
                                   # 'ensemble', 'High_Quality_ChEMBL_compounds',
                                   # 'Small_Molecule_Druggable_Genome_Member',
                                   # 'Bucket_1_ab', 'Bucket_2_ab', 'Bucket_3_ab', 'Bucket_4_ab',
                                   # 'Bucket_5_ab', 'Bucket_6_ab', 'Bucket_7_ab',
                                   # 'Bucket_8_ab', 'Bucket_9_ab', 'Bucket_sum_ab', 'Top_bucket_ab',
                                   # 'Clinical_Precedence_ab', 'Predicted_Tractable_ab_High_confidence', 'Predicted_Tractable_ab_Medium_to_low_confidence', 'Category_ab',
                                   # 'Uniprot_high_conf_loc', 'GO_high_conf_loc',
                                   # 'Uniprot_med_conf_loc',
                                   # 'GO_med_conf_loc', 'Transmembrane', 'Signal_peptide', 'HPA_main_location',
                                   # 'Bucket_1_PROTAC', 'Bucket_2_PROTAC', 'Bucket_3_PROTAC', 'Bucket_4_PROTAC',
                                   # 'Bucket_5_PROTAC', 'Bucket_6_PROTAC', 'Bucket_7_PROTAC',
                                   # 'Bucket_8_PROTAC', 'Bucket_9_PROTAC', 'Bucket_sum_PROTAC', 'Top_bucket_PROTAC',
                                   # 'Bcell_mean', 'NKcell_mean', 'Hepatocytes_mean', 'MouseNeuorons_mean',
                                   # 'Max_halflife',
                                   # 'number_of_ubiquitination_sites',
                                   # 'full_id', 'PROTAC_location_Bucket',
                                   'Bucket_1_othercl', 'Bucket_2_othercl', 'Bucket_3_othercl',
                                   'Bucket_sum_othercl', 'Top_bucket_othercl', 
                                   'drug_chembl_ids_othercl', 'drug_names_othercl', 'clinical_phases_othercl'
                                   ]]

        # self.out_df.sort_values(['Clinical_Precedence', 'Discovery_Precedence', 'Predicted_Tractable'],
        #                         ascending=[False, False, False], inplace=True)

        # Score each category, and label highest category
        self.out_df['Clinical_Precedence_othercl'] = self.out_df.apply(self._clinical_precedence, axis=1)

        self.out_df['Category_othercl'] = 'Unknown'
        self.out_df.loc[(self.out_df['Top_bucket_othercl'] <= 3), 'Category_othercl'] = 'Clinical_Precedence_othercl'
        # self.out_df.loc[(self.out_df['Top_bucket_othercl'] == 4) | (self.out_df['Top_bucket_othercl'] == 5),
        #                 'Category_othercl'] = 'Predicted_Tractable__High_confidence'
        # self.out_df.loc[(self.out_df['Top_bucket_othercl'] == 6) | (self.out_df['Top_bucket_othercl'] == 7) | (
        #     self.out_df['Top_bucket_othercl'] == 8) | (self.out_df['Top_bucket_othercl'] == 9),
        #     'Category_othercl'] = 'Predicted_Tractable__Medium_to_low_confidence'

        # self.out_df = self.out_df[(self.out_df['Top_bucket'] < 9 ) | (self.out_df['Top_bucket_ab'] < 10) | (self.out_df['Top_bucket_othercl'] < 10) ]

        # Cleaning columns
        self.out_df['drug_chembl_ids_othercl'].fillna('', inplace=True)
        self.out_df['drug_names_othercl'].fillna('', inplace=True)
        self.out_df['clinical_phases_othercl'].fillna('', inplace=True)

        # create dictionaries from 'drug_chembl_ids_' and 'clinical_phases_'/'drug_names_'
        self.out_df['drug_names_dict_othercl'] = self.out_df.apply(lambda row : dict(zip(row['drug_chembl_ids_othercl'].split(","), row['drug_names_othercl'].split(" | "))), axis=1)
        self.out_df['clinical_phases_dict_othercl'] = self.out_df.apply(lambda row : dict(zip(row['drug_chembl_ids_othercl'].split(","), row['clinical_phases_othercl'].split(","))), axis=1)

        # Cleaning column: setting selected culumns in list format to improve visualization e.g. with Excel
        # and remove duplicates while keeping order using "list(dict.fromkeys(lst))"
        self.out_df['drug_chembl_ids_othercl'] = self.out_df['drug_chembl_ids_othercl'].apply(lambda x: list(dict.fromkeys(x.split(","))))
        self.out_df['drug_names_othercl'] = self.out_df['drug_names_othercl'].apply(lambda x: list(dict.fromkeys(x.split(" | "))))

        print(self.out_df.columns)

        return self.out_df.astype({x: 'int64' for x in self.out_df.columns if "Bucket" in x})



    @staticmethod
    def othercl2json(d):
        return {
            'Bucket_scores': {'Bucket_1_othercl':d.Bucket_1_othercl, 'Bucket_2_othercl':d.Bucket_2_othercl, 'Bucket_3_othercl':d.Bucket_3_othercl},
            'Bucket_evaluation': {'Bucket_sum_othercl':d.Bucket_sum_othercl, 'Top_bucket_othercl':d.Top_bucket_othercl},
            'Category_scores': {'Clinical_Precedence_othercl':d.Clinical_Precedence_othercl}, 
            'Category_evaluation': {'Top_Category_othercl':d.Category_othercl},
#            'Bucket_evidences': {'Bucket_1_2_3_othercl': {'drug_chembl_ids_othercl':d.drug_chembl_ids_othercl, 'drug_names_othercl':d.drug_names_othercl}}
            'Bucket_evidences': {'Bucket_1_2_3_othercl': {'drugs_in_clinic':d.drug_names_dict_othercl, 'max_clinical_phase':d.clinical_phases_dict_othercl}}
            }

            