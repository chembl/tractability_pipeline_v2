#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  3 10:38:51 2020

@author: Melanie Schneider
"""
# import io
# import re
# import zipfile
# import argparse
# import datetime
import json
import sys
import time
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


from ot_tractability_pipeline_v2.queries_protac import *

DATA_PATH = pkg_resources.resource_filename('ot_tractability_pipeline_v2', 'data/')

class Protac_buckets(object):
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

        # If antibody results are to be combined with small molecule results, append antibody columns to sm results
        # Otherwise, use the id_xref dataframe

        if prev_output is not None:
            go_data = self.id_xref[['ensembl_gene_id', 'go.CC']]
            self.out_df = prev_output.merge(go_data, how='outer', on='ensembl_gene_id')

        else:
            self.out_df = self.id_xref

        # ChEMBL currently not used

        # if database_url is None:
        #     database_url = os.getenv('CHEMBL_DB')
        #
        #
        # # Create ChEMBL DB connection
        # self.engine = create_engine(database_url)


    def _high_conf_locations(self, row):

        if isinstance(row['Uniprot_high_conf_loc'], str):
            if len(row['Uniprot_high_conf_loc']) == 0:
                row['Uniprot_high_conf_loc'] = []
            else: 
                row['Uniprot_high_conf_loc'] = eval(row['Uniprot_high_conf_loc'])

        if isinstance(row['GO_high_conf_loc'], str):
            if len(row['GO_high_conf_loc']) == 0:
                row['GO_high_conf_loc'] = []
            else: 
                row['GO_high_conf_loc'] = eval(row['GO_high_conf_loc'])

        try:
            len(row['Uniprot_high_conf_loc'])
        except TypeError:
            row['Uniprot_high_conf_loc'] = []

        try:
            len(row['GO_high_conf_loc'])
        except TypeError:
            row['GO_high_conf_loc'] = []



        if len(row['Uniprot_high_conf_loc']) == 0 and len(row['GO_high_conf_loc']) == 0 and row['PROTAC_location_Bucket'] == 5:
            return 5




        # locations = [x[0].lower().strip() for x in eval(row['Uniprot_high_conf_loc'])] + [x[0] for x in eval(row['GO_high_conf_loc'])]
        locations = [x[0].lower().strip() for x in row['Uniprot_high_conf_loc']] + [x[0] for x in row['GO_high_conf_loc']]
        accepted_locations = list(set(locations) & set(self.good_locations))
        grey_locations = list(set(locations) & set(self.grey_locations))
        # bad_locations = list(set(locations) & set(self.bad_locations))



        if len(accepted_locations) > 0:
            return 1

        elif len(grey_locations) > 0:
            return 3

        # elif len(bad_loactions) > 0:
        #     return 7


        # If high conf locations are known, but not in self.good_locations or self.grey_locations, they are assumed to
        # be bad

        elif row['PROTAC_location_Bucket'] == 6:
            return 7

        else:
            #print('locations',locations, 'accepted', accepted_locations, 'grey', grey_locations)

            return row['PROTAC_location_Bucket']




    def _med_conf_locations(self,row):

        if isinstance(row['Uniprot_med_conf_loc'], str):
            if len(row['Uniprot_med_conf_loc']) == 0:
                row['Uniprot_med_conf_loc'] = []
            else: 
                row['Uniprot_med_conf_loc'] = eval(row['Uniprot_med_conf_loc'])

        if isinstance(row['GO_med_conf_loc'], str):
            if len(row['GO_med_conf_loc']) == 0:
                row['GO_med_conf_loc'] = []
            else: 
                row['GO_med_conf_loc'] = eval(row['GO_med_conf_loc'])

        try:
            len(row['Uniprot_med_conf_loc'])
        except TypeError:
            row['Uniprot_med_conf_loc'] = []

        try:
            len(row['GO_med_conf_loc'])
        except TypeError:
            row['GO_med_conf_loc'] = []

        # should this be `if len(row['Uniprot_high_conf_loc']) == 0 and len(row['GO_high_conf_loc']) == 0 and row['PROTAC_location_Bucket'] == 5:`?
        if len(row['Uniprot_med_conf_loc']) == 0 and len(row['GO_med_conf_loc']) == 0:
            return 5



        #locations = [x[0].lower().strip() for x in eval(row['Uniprot_med_conf_loc'])] + [x[0] for x in eval(row['GO_med_conf_loc'])]
        locations = [x[0].lower().strip() for x in row['Uniprot_med_conf_loc']] + [x[0] for x in row['GO_med_conf_loc']]
        accepted_locations = list(set(locations) & set(self.good_locations))
        grey_locations = list(set(locations) & set(self.grey_locations))
        # bad_locations = list(set(locations) & set(self.bad_locations))


        if len(accepted_locations) > 0:
            return 2

        elif len(grey_locations) > 0:
            return 4

        # elif len(bad_loactions) > 0:
        #     return 6

        # If high conf locations are known, but not in self.good_locations or self.grey_locations, they are assumed to
        # be bad

        else:
            return 6



    def _PROTAC_location_bucket(self):

        ''''
        For PROTACs, only intracellular targets are suitable. Therefore, we will assign a score based on location to
        allow the targets to be filtered to those in the cytosol, nucleus or membrane with accessible portion.

        1 - High confidence good location
        2 - Med confidence good location
        3 - High confidence grey location
        4 - Med condifence grey location
        5 - Unknown location
        6 - Med confidence bad location
        7 - High confidence bad location
        '''

        self.good_locations = ['cytoplasm', 'cytosol', 'nucleus']
        self.grey_locations = ['membrane']
        # self.bad_locations = ['secreted']



        self.out_df['PROTAC_location_Bucket'] = 0

        self.out_df['PROTAC_location_Bucket'] = self.out_df.apply(self._med_conf_locations, axis=1)
        self.out_df['PROTAC_location_Bucket'] = self.out_df.apply(self._high_conf_locations, axis=1)



    ##############################################################################################################
    #
    # Functions relating to buckets 1-3
    # Clinical PROTAC targets
    #
    ##############################################################################################################

    def _assign_buckets_1_to_3(self):
        '''
        Merge the results of the ChEMBL search with the OT data (right join, to keep all OT targets)

        Group activity data by target, assign the Max Phase for each targets, and use it to assign buckets 1 to 3

        Currently, PROTACs are not labelled in ChEMBL, and at the time of writing, only Androgen Receptor
        (ENSG00000169083) has a phase 1 PROTAC. For now, this will be returned manually, but a decision should be made
        about how PROTACs are labelled in ChEMBL

        :return:
        '''

        self.out_df['Bucket_1_PROTAC'] = 0
        self.out_df['Bucket_2_PROTAC'] = 0
        self.out_df['Bucket_3_PROTAC'] = 0

        self.out_df.loc[(self.out_df['ensembl_gene_id'] == 'ENSG00000169083'), 'Bucket_3_PROTAC'] = 1

    ##############################################################################################################
    #
    # Functions relating to buckets 4 and 5
    # Protein Turnover
    #
    ##############################################################################################################

    def _assign_bucket_4_and_5(self):
        '''
        Protein Turnover
        '''

        self.out_df['Bucket_4_PROTAC'] = 0
        self.out_df['Bucket_5_PROTAC'] = 0

        df = pd.read_csv(os.path.join(DATA_PATH, 'protein_half_life_hq.csv'))

        df = df.merge(self.out_df, right_on='symbol', left_on='gene_name', how='right')
        df = df.groupby('ensembl_gene_id', as_index=False).max(numeric_only=True)
        df['Max_halflife'].fillna(-1, inplace=True)

        self.out_df = df.merge(self.out_df, how='right', on='ensembl_gene_id', suffixes=['_drop', ''])

        self.out_df.loc[(self.out_df['Max_halflife'] >= 24), 'Bucket_4_PROTAC'] = 1
        self.out_df.loc[(self.out_df['Max_halflife'] > 10) & (self.out_df['Max_halflife'] < 24), 'Bucket_5_PROTAC'] = 1

    ##############################################################################################################
    #
    # Functions relating to buckets 6
    # Known ubiquitination sites
    #
    ##############################################################################################################

    def _assign_bucket_6(self):
        '''
        Known ubiquitation sites
        '''

        ub_df = pd.read_csv(os.path.join(DATA_PATH, 'ubiquitination_sites.csv'))
        self.out_df = ub_df.merge(self.out_df, on='symbol', how='right')

        self.out_df['Bucket_6_PROTAC'] = 0
        self.out_df.loc[(self.out_df['number_of_ubiquitination_sites'] > 0), 'Bucket_6_PROTAC'] = 1

    ##############################################################################################################
    #
    # Functions relating to bucket 7
    # Predicted ubiquitination sites
    #
    ##############################################################################################################

    def _assign_bucket_7(self):
        '''
        Predicted ubiquitination sites

        '''
        self.out_df['Bucket_7_PROTAC'] = 0

    ##############################################################################################################
    #
    # Functions relating to buckets 8
    # Taregts mentioned in PROTAC literature
    #
    ##############################################################################################################

    def _search_papers(self):

        url = urllib2.urlopen("https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=%22proteolysis%20targeting%20chimera%22&resultType=lite&cursorMark=*&pageSize=1000&format=json")
        data = url.read()
        try: data = json.loads(data.decode())
        except UnicodeDecodeError: data = json.loads(data)
        df = pd.read_json(json.dumps(data['resultList']['result']), orient='records')

        return df[['authorString', 'id', 'issue',
                   'journalTitle', 'pmcid',
                   'pmid', 'pubType', 'pubYear', 'source', 'title', 'tmAccessionTypeList']]

    def _get_tagged_targets(self):

        articles = list(self.papers_df['search_id'].unique())

        # API only able to accept 8 IDs at a time
        n = 7
        chunks = ['&'.join(articles[i:i + n]) for i in range(0, len(articles), n)]
        df_lists = []
        tags_list = []
        for chunk in chunks:
            url_s = 'https://www.ebi.ac.uk/europepmc/annotations_api/annotationsByArticleIds?{}&type=Gene_Proteins&format=JSON'.format(chunk)
            url = urllib2.urlopen(url_s)
            data = url.read()
            try: data = json.loads(data.decode())
            except UnicodeDecodeError: data = json.loads(data)
            annot_df = json_normalize(data,
                                      record_path='annotations')  # pd.read_json(json.dumps(data), orient='records')
            tags_df = json_normalize(data, record_path=['annotations', 'tags'])
            df_lists.append(annot_df)
            tags_list.append(tags_df)

            time.sleep(1.5)

        self.annotations = pd.concat(df_lists)
        self.tags = pd.concat(tags_list)
        self.annotations.reset_index(inplace=True)
        self.tags.reset_index(inplace=True)

        if self.store_fetched: 
            self.annotations.to_csv("{}/protac_pmc_annotations.csv".format(self.store_fetched), encoding='utf-8')
            self.tags.to_csv("{}/protac_pmc_tags.csv".format(self.store_fetched), encoding='utf-8')

    def _extract_uniprot(self, row):
        try:
            return row['uri'].split('/')[-1]
        except AttributeError:
            return row['uri']

    def _extract_id(self, row):
        try:
            short_id = row['id'].split('/')[-1].split('#')[0]
            return short_id
        except AttributeError:
            return None

    def _process_IDs(self):
        grouped_tags = self.tags.groupby('name').first()
        tagged_annotations = self.annotations.merge(grouped_tags, how='left', left_on='exact', right_on='name')
        tagged_annotations['accession'] = tagged_annotations.apply(self._extract_uniprot, axis=1)

        tagged_annotations['short_id'] = tagged_annotations.apply(self._extract_id, axis=1)

        # tagged_annotations['short_id']
        joined = tagged_annotations.merge(self.papers_df, left_on='short_id', right_on='id', how='inner')

        return joined[['accession', 'prefix', 'exact', 'postfix', 'section', 'full_id', 'journalTitle']]

    def _search_ID(self, row):
        return "articleIds={}%3A{}".format(row['source'], row['id'])

    def _full_ID(self, row):
        return "http://europepmc.org/abstract/{}/{}#eur...".format(row['source'], row['id'])

    def _assign_bucket_8(self):
        '''
        Mentioned in PROTAC literature
        '''

        self.papers_df = self._search_papers()
        
        if self.store_fetched: 
            self.papers_df.to_csv("{}/protac_pmc_papers.csv".format(self.store_fetched), encoding='utf-8')

        self.papers_df['search_id'] = self.papers_df.apply(self._search_ID, axis=1)
        self.papers_df['full_id'] = self.papers_df.apply(self._full_ID, axis=1)

        self._get_tagged_targets()

        tagged_targets_df = self._process_IDs()

        self.out_df = self.out_df.merge(tagged_targets_df, how='left', on='accession')

        self.out_df['Bucket_8_PROTAC'] = 0
        self.out_df.loc[(~self.out_df['full_id'].isna()), 'Bucket_8_PROTAC'] = 1

    ##############################################################################################################
    #
    # Functions relating to buckets 8
    # Small Molecule Tractable
    #
    ##############################################################################################################

    def _assign_bucket_9(self):
        '''
        Small molecule tractable
        '''

        self.out_df['Bucket_9_PROTAC'] = 0
        self.out_df.loc[(self.out_df['Top_bucket'] < 9), 'Bucket_9_PROTAC'] = 1

    ##############################################################################################################
    #
    # Higher level functions relating to the overall process
    #
    #
    ##############################################################################################################

    # def _clinical_precedence(self, s):
    #     return 1 * s['Bucket_1_ab'] + 0.7 * s['Bucket_2_ab'] + 0.2 * s['Bucket_3_ab']
    #
    # def _high_conf_pred(self, s):
    #     return 0.7 * s['Bucket_4_ab'] + 0.3 * s['Bucket_5_ab']
    #
    # def _med_conf_pred(self, s):
    #     return 0.4 * s['Bucket_6_ab'] + 0.25 * s['Bucket_7_ab'] + 0.25 * s['Bucket_8_ab'] + 0.1 * s['Bucket_9_ab']

    def _summarise_buckets(self):

        self.out_df['Top_bucket_PROTAC'] = 10
        for x in range(9, 0, -1):
            self.out_df.loc[(self.out_df['Bucket_{}_PROTAC'.format(x)] == 1), 'Top_bucket_PROTAC'] = x
            self.out_df['Bucket_{}_PROTAC'.format(x)].fillna(0, inplace=True)

        self.out_df['Bucket_sum_PROTAC'] = self.out_df['Bucket_1_PROTAC'] + self.out_df['Bucket_2_PROTAC'] + \
                                           self.out_df[
                                               'Bucket_3_PROTAC'] + self.out_df['Bucket_4_PROTAC'] + self.out_df[
                                               'Bucket_5_PROTAC'] + self.out_df['Bucket_6_PROTAC'
                                           ] + self.out_df['Bucket_7_PROTAC'] + self.out_df['Bucket_8_PROTAC'] + \
                                           self.out_df[
                                               'Bucket_9_PROTAC']

        self.out_df.set_index('ensembl_gene_id')

    def assign_buckets(self):
        '''
        Assigns the supplied list of gene IDs into their corresponding tractability buckets.
        :return: A Pandas DataFrame containing the Ensembl gene ID and associated tractability bucket
        '''

        self._PROTAC_location_bucket()
        self._assign_buckets_1_to_3()
        self._assign_bucket_4_and_5()
        self._assign_bucket_6()
        self._assign_bucket_7()
        self._assign_bucket_8()
        self._assign_bucket_9()

        self._summarise_buckets()

        # try:

        # self.out_df.index = self.out_df['ensembl_gene_id']
        self.out_df = self.out_df.groupby('ensembl_gene_id').first()

        # Columns to keep. This includes columns from the small molecule pipeline
        self.out_df = self.out_df[['accession', 'symbol',
                                   'Bucket_1', 'Bucket_2', 'Bucket_3', 'Bucket_4',
                                   'Bucket_5', 'Bucket_6', 'Bucket_7',
                                   'Bucket_8', 'Bucket_sum', 'Top_bucket', 'Category',
                                   'Clinical_Precedence', 'Discovery_Precedence', 'Predicted_Tractable',
                                   'PDB_Known_Ligand',
                                   'ensemble', 'High_Quality_ChEMBL_compounds',
                                   'Small_Molecule_Druggable_Genome_Member',
                                   'Bucket_1_ab', 'Bucket_2_ab', 'Bucket_3_ab', 'Bucket_4_ab',
                                   'Bucket_5_ab', 'Bucket_6_ab', 'Bucket_7_ab',
                                   'Bucket_8_ab', 'Bucket_9_ab', 'Bucket_sum_ab', 'Top_bucket_ab',
                                   'Clinical_Precedence_ab', 'Predicted_Tractable__High_confidence', 'Predicted_Tractable__Medium_to_low_confidence', 'Category_ab',
                                   'Uniprot_high_conf_loc', 'GO_high_conf_loc',
                                   'Uniprot_med_conf_loc',
                                   'GO_med_conf_loc', 'Transmembrane', 'Signal_peptide', 'HPA_main_location',
                                   'Bucket_1_PROTAC', 'Bucket_2_PROTAC', 'Bucket_3_PROTAC', 'Bucket_4_PROTAC',
                                   'Bucket_5_PROTAC', 'Bucket_6_PROTAC', 'Bucket_7_PROTAC',
                                   'Bucket_8_PROTAC', 'Bucket_9_PROTAC', 'Bucket_sum_PROTAC', 'Top_bucket_PROTAC',
                                   'Bcell_mean', 'NKcell_mean', 'Hepatocytes_mean', 'MouseNeuorons_mean',
                                   'Max_halflife',
                                   'number_of_ubiquitination_sites',
                                   'full_id', 'PROTAC_location_Bucket'
                                   ]]
        self.out_df.sort_values(['Clinical_Precedence', 'Discovery_Precedence', 'Predicted_Tractable'],
                                ascending=[False, False, False], inplace=True)


        self.out_df = self.out_df[(self.out_df['Top_bucket'] < 9) | (self.out_df['Top_bucket_ab'] < 10) | (
                    self.out_df['Top_bucket_PROTAC'] < 10)]

        return self.out_df.astype({x: 'int64' for x in self.out_df.columns if "Bucket" in x})
