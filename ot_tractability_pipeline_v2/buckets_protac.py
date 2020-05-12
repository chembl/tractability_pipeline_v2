#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script contains all functions for the PROTAC workflow.

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
import numpy as np
import pandas as pd
import pkg_resources
import ast
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

        # ChEMBL DB connection
        self.engine = Pipeline_setup.engine

        # # use function from Pipeline setup
        # self.make_request = Pipeline_setup.make_request

        # All chembl data loaded into here
        self.all_chembl_targets = None
        # Human proteome data loaded into here
        self.human_proteome = None

        # If antibody results are to be combined with small molecule results, append antibody columns to sm results
        # Otherwise, use the id_xref dataframe

        if prev_output is not None:
            go_data = self.id_xref[['ensembl_gene_id', 'go.CC']]
            self.out_df = prev_output.merge(go_data, how='outer', on='ensembl_gene_id')

        else:
            self.out_df = self.id_xref
                        
            '''
            Uniprot (loc)
            '''
            try:
                # try reading the output from store_tetched antibody pipeline to avoid getting data again
                uniprot_loc_data = pd.read_csv("{}/ab_uniprot_locations_processed.csv".format(self.store_fetched))
                # self.out_uniprot_loc_data = self.out_uniprot_loc_data.merge(uniprot_loc_data, how='left', on='accession')

            except IOError:
                print("--- \nWarning in PROTAC workflow: Uniprot location file 'ab_uniprot_locations_processed.csv' is not available \
                      from antibody workflow. (see section 'Uniprot (loc)' in buckets_ab.py) \n-> trying to reasses from source \nThis will take some time ...\n---")
                uniprot_loc_data = None
            
            # =============================================================================
            #     The following code is only executed if no previously calculated ab_uniprot_locations_processed.csv is found     
            # =============================================================================
            if uniprot_loc_data is None:
                
                from ot_tractability_pipeline_v2.buckets_ab import Antibody_buckets
                # make_request = Antibody_buckets.make_request
                # post_request_uniprot = Antibody_buckets.post_request_uniprot
                # split_loc = Antibody_buckets.split_loc
                # _check_evidence = Antibody_buckets._check_evidence

                            
                # Return all reviewed and Human targets
                url = "uniprot/?format=tab&query=*&fil=reviewed%3ayes+AND+organism%3a%22Homo+sapiens+(Human)+%5b9606%5d%22&columns=id,comment(SUBCELLULAR+LOCATION),comment(DOMAIN),feature(DOMAIN+EXTENT),feature(INTRAMEMBRANE),feature(TOPOLOGICAL+DOMAIN),feature(TRANSMEMBRANE),feature(SIGNAL)"
                data = ['P42336', 'P60484']
                location = Antibody_buckets.post_request_uniprot(url, data)
                location = [x.split('\t') for x in location.split('\n')]
                uniprot_loc_data = pd.DataFrame(location[1:], columns=location[0])
                uniprot_loc_data['uniprot_loc_test'] = uniprot_loc_data['Subcellular location [CC]']
                uniprot_loc_data['Subcellular location [CC]'] = uniprot_loc_data['Subcellular location [CC]'].apply(Antibody_buckets.split_loc)
        
                uniprot_loc_data.rename(columns={'Entry': 'accession'}, inplace=True)
        
                def get_uniprot_high_conf(s):            
                    all_uniprot_high_conf = [(a[1], a[0]) for a in s['Subcellular location [CC]'] if Antibody_buckets._check_evidence(a[0])]    
                    return all_uniprot_high_conf

                def get_uniprot_med_conf(s):            
                    all_uniprot_med_conf = [(a[1], a[0]) for a in s['Subcellular location [CC]'] if not Antibody_buckets._check_evidence(a[0])]
                    return all_uniprot_med_conf
        
                uniprot_loc_data['Uniprot_high_conf_loc'] = zip(*uniprot_loc_data.apply(get_uniprot_high_conf, axis=1))
                uniprot_loc_data['Uniprot_med_conf_loc'] = zip(*uniprot_loc_data.apply(get_uniprot_med_conf, axis=1))
                # uniprot_loc_data['Uniprot_high_conf_loc'] = zip(*uniprot_loc_data.apply(Antibody_buckets.get_uniprot_high_conf, axis=1))
                # uniprot_loc_data['Uniprot_med_conf_loc'] = zip(*uniprot_loc_data.apply(Antibody_buckets.get_uniprot_med_conf, axis=1))
    
                if self.store_fetched: 
                    uniprot_loc_data.to_csv("{}/protac_uniprot_location.csv".format(self.store_fetched))

            else:
                print("\t--- Uniprot location data used from Antibody workflow ---")

            # finally:
            print(uniprot_loc_data.columns)
            self.out_df = self.out_df.merge(uniprot_loc_data, how='left', on='accession')
            # print(self.out_df.columns)

            '''
            GO CC (loc)
            '''
            try:
                # try reading the output from store_tetched antibody pipeline to avoid getting data again
                go_loc_data = pd.read_csv("{}/ab_GO_locations_processed.csv".format(self.store_fetched))
                print("\t--- GO location data used from Antibody workflow file ab_GO_locations_processed.csv ---")
                self.out_df = self.out_df.merge(go_loc_data, how='left', on='accession')
                
            except IOError:
                print("--- \nWarning in PROTAC workflow: GO location file 'ab_GO_locations_processed.csv' is not available \
                      from antibody workflow. (see section 'GO CC' in buckets_ab.py) \n-> trying to reasses from source \nThis will take some time ...\n---")
                go_loc_data = None
       
            # =============================================================================
            #     The following code is only executed if no previously calculated ab_GO_locations_processed.csv is found     
            # =============================================================================
            if go_loc_data is None:
                # go_loc_data = self.out_df
                
                def get_go_loc(s):
                    try:
                        cc = s['go.CC']
                    except:
                        # return 0, [], 0, []
                        return [], []
            
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
                        # return 0, [], 0, []
                        return [], []
            
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
            
                        if go_id in Antibody_buckets.accepted_go_locs.keys():
                            if confidence == 'High':
                                accepted_high_conf_loc.append(go_loc)
                            elif confidence == 'Medium':
                                accepted_med_conf_loc.append(go_loc)
            
                    # b5_flag = 0
                    # b8_flag = 0
            
                    # if len(accepted_high_conf_loc) > 0:
                    #     b5_flag = 1
                    # elif len(accepted_med_conf_loc) > 0:
                    #     b8_flag = 1
            
                    # return b5_flag, high_conf_loc, b8_flag, med_conf_loc
                    return high_conf_loc, med_conf_loc
            
                '''
                GO CC
                '''
                # self.out_df['Bucket_5_ab'], self.out_df['GO_high_conf_loc'], self.out_df['Bucket_8_ab'], self.out_df[
                #     'GO_med_conf_loc'] = zip(*self.out_df.apply(self.get_go_loc, axis=1))
                # go_loc_data['GO_high_conf_loc'], go_loc_data['GO_med_conf_loc'] = zip(*go_loc_data.apply(get_go_loc, axis=1))
                self.out_df['GO_high_conf_loc'], self.out_df['GO_med_conf_loc'] = zip(*self.out_df.apply(get_go_loc, axis=1))
                
                # save processed GO location info to file
                if self.store_fetched: 
                    # go_loc_data.loc[:, go_loc_data.columns.isin(['accession', 'GO_high_conf_loc', 'GO_med_conf_loc'])].to_csv(
                    #     "{}/protac_GO_locations_processed.csv".format(self.store_fetched))
                    self.out_df.loc[:, self.out_df.columns.isin(['accession', 'GO_high_conf_loc', 'GO_med_conf_loc'])].to_csv(
                        "{}/protac_GO_locations_processed.csv".format(self.store_fetched))
        
                # print(go_loc_data.columns)
                print(self.out_df.columns)
                

            else:
                # print("\t--- GO location data used from Antibody workflow ---")
                print(go_loc_data.columns)
                print(self.out_df.columns)



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

        # If high conf locations are known, but not in self.good_locations or self.grey_locations, they are assumed to be bad
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

        # If high conf locations are known, but not in self.good_locations or self.grey_locations, they are assumed to be bad
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

        print("\t- Assessing PROTAC location bucket...")

        self.good_locations = ['cytoplasm', 'cytosol', 'nucleus']
        self.grey_locations = ['membrane']
        # self.bad_locations = ['secreted']


        self.out_df['PROTAC_location_Bucket'] = 0

        self.out_df['PROTAC_location_Bucket'] = self.out_df.apply(self._med_conf_locations, axis=1)
        self.out_df['PROTAC_location_Bucket'] = self.out_df.apply(self._high_conf_locations, axis=1)

        print(self.out_df.columns)


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
        print("\t- Assessing clinical buckets 1-3...")

        self.out_df['Bucket_1_PROTAC'] = 0
        self.out_df['Bucket_2_PROTAC'] = 0
        self.out_df['Bucket_3_PROTAC'] = 0

        # AR PROTAC ARV-110 is currently in clinical trials for metastatic castration-resistant prostate cancer. 
        self.out_df.loc[(self.out_df['ensembl_gene_id'] == 'ENSG00000169083'), 'Bucket_3_PROTAC'] = 1
        # ER PROTAC ARV-471 is currently in clinical trials for advanced or metastatic ER+/Her2− breast cancer. 
        self.out_df.loc[(self.out_df['ensembl_gene_id'] == 'ENSG00000091831'), 'Bucket_3_PROTAC'] = 1

        # Cleaning column: setting selected culumns in list format to improve visualization e.g. with Excel
        # and remove duplicates while keeping order using "list(dict.fromkeys(lst))"
        # self.out_df['drug_chembl_ids_PROTAC'].fillna('', inplace=True)
        # self.out_df['drug_chembl_ids_PROTAC'] = self.out_df['drug_chembl_ids_PROTAC'].apply(lambda x: list(dict.fromkeys(x.split(","))))

        print(self.out_df.columns)


    ##############################################################################################################
    #
    # Functions relating to bucket 4
    # Protein Turnover
    #
    ##############################################################################################################

    # def _assign_bucket_4_and_5(self):
    def _assign_bucket_4(self):
        '''
        Protein Turnover
        '''

        print("\t- Assessing protein turnover bucket 4...")

        self.out_df['Bucket_4_PROTAC'] = 0
        # self.out_df['Bucket_5_PROTAC'] = 0

        df = pd.read_csv(os.path.join(DATA_PATH, 'protein_half_life_hq.csv'))

        df = df.merge(self.out_df, right_on='symbol', left_on='gene_name', how='right')
        df = df.groupby('ensembl_gene_id', as_index=False).max(numeric_only=True)
        df['Max_halflife'].fillna(-1, inplace=True)
        df['Min_halflife'].fillna(-1, inplace=True)

        self.out_df = df.merge(self.out_df, how='right', on='ensembl_gene_id', suffixes=['_drop', ''])
        
        # Assign bucket if Half-life data is available
        self.out_df.loc[(self.out_df['Min_halflife'] > 0), 'Bucket_4_PROTAC'] = 1
        # self.out_df.loc[(self.out_df['Max_halflife'] >= 24), 'Bucket_4_PROTAC'] = 1
        # self.out_df.loc[(self.out_df['Max_halflife'] > 10) & (self.out_df['Max_halflife'] < 24), 'Bucket_5_PROTAC'] = 1

        print(self.out_df.columns)


    ##############################################################################################################
    #
    # Functions relating to bucket 5
    # UniProt ubiquitination indication based on keyword:"Ubl conjugation [KW-0832]" 
    #
    ##############################################################################################################

    #staticmethod
    # def get_uniprot_ubl(self):
    def _assign_bucket_5(self):
        '''
        UniProt ubiquitination indication based on keyword:"Ubl conjugation [KW-0832]"
        '''

        print("\t- Assessing UniProt ubiquitination indication bucket 5...")

        full_url = 'https://www.uniprot.org/uniprot/?query=keyword%3A%22Ubl+conjugation+%5BKW-0832%5D%22+organism%3A%22Homo+sapiens+%28Human%29+%5B9606%5D%22&format=tab&columns=id,comment(PTM),feature(CROSS%20LINK)' #,database(PhosphoSitePlus)
        
        from ot_tractability_pipeline_v2.buckets_ab import Antibody_buckets
        Uniprot_ubl_conjugation = Antibody_buckets.make_request(full_url, data=None)
        # url = 'uniprot/?query=keyword%3A%22Ubl+conjugation+%5BKW-0832%5D%22+organism%3A%22Homo+sapiens+%28Human%29+%5B9606%5D%22&format=tab&columns=id,comment(PTM),database(PhosphoSitePlus)'
        # data = ['P42336', 'P60484']
        # Uniprot_ubl_conjugation = Antibody_buckets.post_request_uniprot(url, data)

        Uniprot_ubl_conjugation = [x.split('\t') for x in Uniprot_ubl_conjugation.split('\n')]
        df = pd.DataFrame(Uniprot_ubl_conjugation[1:], columns=Uniprot_ubl_conjugation[0])

        df.rename(columns={'Entry': 'accession', 
                           'Post-translational modification': 'Uniprot_PTM',
                           'Cross-link': 'Uniprot_CrossLink'
                           }, inplace=True)


        if self.store_fetched: 
            df.to_csv("{}/protac_uniprot_ubl_conjugation.csv".format(self.store_fetched))

        df['Bucket_5_PROTAC'] = 1
        df['Uniprot_keyword'] = 'Yes'

        self.out_df = self.out_df.merge(df, how='left', on='accession')
        # replace NaN with 0
        self.out_df['Bucket_5_PROTAC'].fillna(0, inplace=True)
        self.out_df['Uniprot_keyword'].fillna('No', inplace=True)

        print(self.out_df.columns)

    

    ##############################################################################################################
    #
    # Functions relating to bucket 6
    # Reported ubiquitination sites in PhosphoSitePlus, mUbiSiDa and [Woong et al. 2011]
    #
    ##############################################################################################################

    def _assign_bucket_6(self):
        '''
        Reported ubiquitation sites in PhosphoSitePlus, mUbiSiDa and [Woong et al. 2011]
        '''

        print("\t- Assessing ubiquitation site bucket 6...")
        
        PhosphoSitePlus = pd.read_csv(os.path.join(DATA_PATH, "PSP_Ubiquitination_site_dataset"), skiprows=3, sep='\t')
        mUbiSiDa = pd.read_excel(os.path.join(DATA_PATH, "mUbiSiDa_data_2013_10_22.xlsx"), header=None)
        ub_df = pd.read_csv(os.path.join(DATA_PATH, 'ubiquitination_sites.csv'))
        self.out_df = ub_df.merge(self.out_df, on='symbol', how='right')

        self.out_df['Ub_PhosphoSitePlus'] = 0
        self.out_df['Ub_mUbiSiDa_2013'] = 0
        self.out_df['Bucket_6_PROTAC'] = 0
        self.out_df.loc[(self.out_df['accession'].isin(PhosphoSitePlus['ACC_ID'])), 'Ub_PhosphoSitePlus'] = 1
        self.out_df.loc[(self.out_df['accession'].isin(PhosphoSitePlus['ACC_ID'])), 'Bucket_6_PROTAC'] = 1
        self.out_df.loc[(self.out_df['accession'].isin(mUbiSiDa[0])), 'Ub_mUbiSiDa_2013'] = 1
        self.out_df.loc[(self.out_df['accession'].isin(mUbiSiDa[0])), 'Bucket_6_PROTAC'] = 1
        self.out_df.loc[(self.out_df['number_of_ubiquitination_sites'] > 0), 'Bucket_6_PROTAC'] = 1

        print(self.out_df.columns)



    ##############################################################################################################
    #
    # Functions relating to bucket 7
    # Taregts mentioned in PROTAC literature
    #
    ##############################################################################################################

    def _search_papers(self):
        
#        url = urllib2.urlopen("https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=%22proteolysis%20targeting%20chimera%22&resultType=lite&cursorMark=*&pageSize=1000&format=json")
        # another search query: (("proteolysis targeting chimera") AND ("PROTAC")) OR (("PROTAC") AND (("degradation") OR ("degrade")))
        # where "PROTAC" is present in both options to avoid papers mentioning just once "proteolysis targeting chimera" (e.g. as technique or example)

        # keyword search only in Abstract (keywords in "" need additional 'ABSTRACT:'): 
        # ABSTRACT:((ABSTRACT:"proteolysis targeting chimera" OR ABSTRACT:"proteolysis targeting chimeric") OR (PROTAC AND (degradation OR degrade OR ubiquitin OR proteolysis)))
        url = urllib2.urlopen("https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=ABSTRACT%3A%28%28ABSTRACT%3A%22proteolysis%20targeting%20chimera%22%20OR%20ABSTRACT%3A%22proteolysis%20targeting%20chimeric%22%29%20OR%20%28PROTAC%20AND%20%28degradation%20OR%20degrade%20OR%20ubiquitin%20OR%20proteolysis%29%29%29&resultType=lite&cursorMark=*&pageSize=1000&format=json")
        
        data = url.read()
        try: data = json.loads(data.decode())
        except UnicodeDecodeError: data = json.loads(data)
        df = pd.read_json(json.dumps(data['resultList']['result']), orient='records')
        
        # remove entries with missing 'authorString'
        df = df.loc[~df['authorString'].isna()]
        
        return df[['authorString', 'id', 'issue',
                   'journalTitle', 'pmcid',
                   'pmid', 'pubType', 'pubYear', 'source', 'title', 'tmAccessionTypeList']]

    def _search_ID(self, row):
        return "articleIds={}%3A{}".format(row['source'], row['id'])

    def _full_ID(self, row):
        return "http://europepmc.org/abstract/{}/{}#eur...".format(row['source'], row['id'])

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
            annot_df = pd.json_normalize(data, record_path='annotations')  # pd.read_json(json.dumps(data), orient='records')
            tags_df = pd.json_normalize(data, record_path=['annotations', 'tags'])
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
        

    @staticmethod
    def _extract_uniprot(row):
        try:
            return row['uri'].split('/')[-1]
        except AttributeError:
            return row['uri']

    @staticmethod
    def _extract_id(row):
        try:
            short_id = row['id'].split('/')[-1].split('#')[0]
            return short_id
        except AttributeError:
            return None

    def _process_IDs(self):
#        grouped_tags = self.tags.groupby('name').first()
#        tagged_annotations = self.annotations.merge(grouped_tags, how='left', left_on='exact', right_on='name')
        tagged_annotations = self.annotations.merge(self.tags, how='left', left_index=True, right_index=True)

        # keep only tagged_annotations from Abstract of publications
        tagged_annotations = tagged_annotations.loc[tagged_annotations['section'].astype(str).str.contains('Abstract')]
        
        tagged_annotations['EuropePMC_accession'] = tagged_annotations.apply(self._extract_uniprot, axis=1)
        
        tagged_annotations['short_id'] = tagged_annotations.apply(self._extract_id, axis=1)
        
        tagged_annotations.rename(columns={'id':'link_id'}, inplace=True)
        
        # join tagged_annotations with self.papers_df
        joined = tagged_annotations.merge(self.papers_df, left_on='short_id', right_on='id', how='inner')
        
        # Getting the human proteome from UniProt
        full_url = 'https://www.uniprot.org/uniprot/?query=proteome:UP000005640&format=tab&columns=id,entry%20name,protein%20names,genes'
        from ot_tractability_pipeline_v2.buckets_ab import Antibody_buckets
        Uniprot_human_proteome = Antibody_buckets.make_request(full_url, data=None)
        Uniprot_human_proteome = [x.split('\t') for x in Uniprot_human_proteome.split('\n')]
        human_proteome = pd.DataFrame(Uniprot_human_proteome[1:], columns=Uniprot_human_proteome[0])
        human_proteome.rename(columns={'Entry': 'accession'}, inplace=True)
        # only keep row when 'Entry name' is available (discard NAN row)
        human_proteome = human_proteome.loc[human_proteome['Entry name'].notna()]
        # create 'symbol' column
        human_proteome[['symbol','Human']] = human_proteome['Entry name'].str.split("_",expand=True)
        # create 'gene_name' column (using first entry in 'Gene names')
        human_proteome['gene_name'] = human_proteome['Gene names'].str.split(" ",expand=True)[0]
        # create 'protein_name' column (using primary entry, before names in parentheses in 'Protein names', escaping with \\ is reqired)
        human_proteome['protein_name'] = human_proteome['Protein names'].str.split(" \\(",expand=True)[0]
        
        # as all protein isoforms have different UniProt IDs, only the first occurence of gene_name is kept 
        # (which should be the primary UniProtID) count: 20487
        human_proteome.drop_duplicates(subset="gene_name", keep='first', inplace=True)

        # Join on human_proteome data to get human UniProt accession IDs:
        # 1. join with human_proteome on 'symbol' (higher chances of getting the correct/dominant isoform)
        # 2. join with human_proteome on 'gene_name' (occurs several times, for all isoforms)
        # 3. join both new datasets
        tagged_targets_on_symbol = joined.merge(human_proteome, left_on='name', right_on='symbol', how='left')
        tagged_targets_on_gene_name = joined.merge(human_proteome, left_on='name', right_on='gene_name', how='left')
        tagged_targets_df = tagged_targets_on_symbol.merge(tagged_targets_on_gene_name[['name','accession','gene_name']], on='name', how='inner')
        # 4. convert name and protein_name to lowercase for merging 
        #(to avoid missing through different uppercase usage in text)
        # 5. join with human_proteome on 'protein_name_lower' (to get full name annotations) = third new dataset 'tagged_targets_on_protein_name'
        # 6. join tagged_targets_df with third new dataset 'tagged_targets_on_protein_name' using lowercase names
        joined['name_lower'] = joined['name'].str.lower()
        human_proteome['protein_name_lower'] = human_proteome['protein_name'].str.lower()
        tagged_targets_on_protein_name = joined.merge(human_proteome, left_on='name_lower', right_on='protein_name_lower', how='left')
        tagged_targets_df = tagged_targets_df.merge(tagged_targets_on_protein_name[['name','accession','protein_name']], on='name', how='left')

        # 7. remove duplicated rows (keep='first' is default), as lists are contained in df:
        # convert the df to str type astype(str), drop duplicates and then select the rows from original df, thus output df still contains lists
        tagged_targets_df = tagged_targets_df.loc[tagged_targets_df.astype(str).drop_duplicates().index]
        
        # rename last joined accession column to "accession_z"
        tagged_targets_df.rename(columns={"accession": "accession_z"}, inplace=True)
        
        # add additional mapping 'accession_z2' column based on protein name from publication matching end of UniProts full protein name
        tagged_targets_df['accession_z2'] = ''
        tagged_targets_df['accession_z2'] = tagged_targets_df['name'].str.lower().apply(lambda x: human_proteome[human_proteome['protein_name_lower'].str.endswith(x)]['accession'].any(0))
        tagged_targets_df['accession_z2'].replace(to_replace=False, value=np.nan, inplace=True)
        
        
        # as new 'accession' column take 'accession_x' (from tagged_targets_on_symbol)
        # if 'accession_x' (from tagged_targets_on_symbol) is NA, take 'accession_y' (from tagged_targets_on_gene_name)
        tagged_targets_df['accession'] = tagged_targets_df['accession_x'].fillna(value=tagged_targets_df['accession_y'])
        # if 'accession' is still NA, take 'accession_z' (from tagged_targets_on_protein_name)
        tagged_targets_df['accession'] = tagged_targets_df['accession'].fillna(value=tagged_targets_df['accession_z'])
        # if 'accession' is still NA, take 'accession_z2' (from 'protein_name_lower'.str.endswith('exact_lower') = substring search at end of full UniProt protein_name_lower)
        tagged_targets_df['accession'] = tagged_targets_df['accession'].fillna(value=tagged_targets_df['accession_z2'])
        # if 'accession' is still NA, take 'EuropePMC_accession'
        tagged_targets_df['accession'] = tagged_targets_df['accession'].fillna(value=tagged_targets_df['EuropePMC_accession'])

        if self.store_fetched: 
            tagged_targets_df.to_csv("{}/protac_pmc_tagged_targets.csv".format(self.store_fetched), encoding='utf-8')

        self.human_proteome = human_proteome[['accession','symbol','gene_name','protein_name']]
        
#        return joined[['accession', 'prefix', 'exact', 'postfix', 'section', 'full_id', 'journalTitle', 'title']]
        return tagged_targets_df[['accession', 'name', 'short_id', 'full_id', 'title']]


    def _assign_bucket_7(self):
        '''
        Mentioned in PROTAC literature
        '''

        print("\t- Assessing PROTAC literature bucket 7...")

        self.papers_df = self._search_papers()
        
        self.papers_df['search_id'] = self.papers_df.apply(self._search_ID, axis=1)
        self.papers_df['full_id'] = self.papers_df.apply(self._full_ID, axis=1)

        if self.store_fetched: 
            self.papers_df.to_csv("{}/protac_pmc_papers.csv".format(self.store_fetched), encoding='utf-8')

        self._get_tagged_targets()

        tagged_targets_df = self._process_IDs()

        def set_strings(x):
            ''' concatenate in string and include only if it is a string (not nan), and exists '''
            return ",".join([y for y in x if isinstance(y,str) and y])
        def set_title_strings(x):
            ''' concatenate in string and include only if it is a string (not nan), and exists '''
            return ";".join([y for y in x if isinstance(y,str) and y])
        
        f = {x: 'first' for x in tagged_targets_df.columns}
        f['short_id'] = set_strings
        f['full_id'] = set_strings
        f['title'] = set_title_strings
        
        tagged_targets_df_grouped = tagged_targets_df.groupby(['accession']).agg(f).reset_index(drop=True)
        
        # removing duplicates        
        def clean_column(df, col):
            df[col] = df[col].apply(lambda x: list(dict.fromkeys(x.split(","))))
        def clean_title_column(df, col):
            df[col] = df[col].apply(lambda x: list(dict.fromkeys(x.split(";"))))
        
        clean_column(tagged_targets_df_grouped, 'short_id')
        clean_column(tagged_targets_df_grouped, 'full_id')
        clean_title_column(tagged_targets_df_grouped, 'title')
        

        self.out_df = self.out_df.merge(tagged_targets_df_grouped, how='left', on='accession')
        
        self.out_df['Bucket_7_PROTAC'] = 0
        self.out_df.loc[(~self.out_df['full_id'].isna()), 'Bucket_7_PROTAC'] = 1
        
        # count literature per target and append as column
        self.out_df['literature_count_PROTAC'] = [len(x) for x in ast.literal_eval(",".join(self.out_df['full_id'].fillna("['']",inplace=False).astype(str).replace(to_replace="['']",value="''",inplace=False)))]

        # adding protein_name from human_proteome to self.out_df
        self.out_df = self.out_df.merge(self.human_proteome[['accession','protein_name']], how='outer', on='accession')        
        # for debugging:
        # self.out_df.to_csv("{}/protac_out_df_checkpoint.csv".format(self.store_fetched), encoding='utf-8')
        print(self.out_df.columns)

    ##############################################################################################################
    #
    # Functions relating to bucket 8
    # Existing active molecule
    # (before: Small Molecule Tractable)
    #
    ##############################################################################################################
    # def _assign_bucket_9(self):
    #     '''
    #     Small molecule tractable
    #     '''
    #     self.out_df['Bucket_9_PROTAC'] = 0
    #     self.out_df.loc[(self.out_df['Top_bucket_sm'] < 9), 'Bucket_9_PROTAC'] = 1


    def _search_chembl_smallmol_active(self):
        '''
        Search for all targets in ChEMBL, get all small molecules and their targets with a pChEMBL >=5
        (see queries_protac.py)
        :return:
        '''

        # print("\t\t- Querying ChEMBL...")
        # small_mol_info = pd.read_sql_query(chembl_small_mol, self.engine)
        self.all_chembl_targets = pd.read_sql_query(chembl_small_mol_active_targets, self.engine)
        # self.all_chembl_targets = self.all_chembl_targets.merge(small_mol_info, on='parent_molregno')

        if self.store_fetched: 
            self.all_chembl_targets.to_csv("{}/protac_all_chembl_active_targets.csv".format(self.store_fetched))

    def _process_protein_complexes(self):
        '''
        For protein complexes, see if we know the binding subunit, and only keep these
        :return:
        '''

        # print("\t\t- Processing protein complexes...")
        pc = self.all_chembl_targets[self.all_chembl_targets['target_type'].str.contains("PROTEIN COMPLEX")]
        not_pc = self.all_chembl_targets[~self.all_chembl_targets['target_type'].str.contains("PROTEIN COMPLEX")]

        n = 1000
        targets = pc['tid'].unique()
        chunks = [targets[i:i + n] for i in range(0, len(targets), n)]

        df_list = []

        # Check if binding sites are defined
        for chunk in chunks:
            q = '''
            select distinct bs.site_id, td.tid 
            from {0}.target_dictionary td, {0}.binding_sites bs
            where td.tid = bs.tid and td.tid IN {1}
            '''.format(CHEMBL_VERSION, tuple(chunk))
            df_list.append(pd.read_sql_query(q, self.engine))

        # Merge will set those with unknown binding site as NAN
        binding_site_info = pd.concat(df_list, sort=False)

        if self.store_fetched: 
            binding_site_info.to_csv("{}/protac_smallmol_chembl_binding_site_info.csv".format(self.store_fetched))

        pc = pc.merge(binding_site_info, how='left', on='tid')
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
            df_list2.append(pd.read_sql_query(q2, self.engine))

        binding_subunit = pd.concat(df_list2, sort=False)

        if self.store_fetched: 
            binding_subunit.to_csv("{}/protac_smallmol_chembl_binding_subunit.csv".format(self.store_fetched))

        temp_pc = pc.merge(binding_subunit, on='accession')
        binding_subunits = temp_pc[temp_pc['component_id'].notnull()]

        self.all_chembl_targets = pd.concat([binding_subunits, undefined, not_pc], sort=False)

    def _assign_bucket_8(self):
        '''
        Target with SM in ChEMBL with activity on target based assay
        Merge the results of the ChEMBL search with the OT data (right join, to keep all OT targets)
        Group activity data by target, add 'compound_chembl_id' and 'pchembl_value' as tupels to target,
        and assign bucket 8
        assign the Max Phase for each targets, and use it to assign buckets 1 to 3
        '''

        print("\t- Assessing existing active molecule bucket 8...")

        def set_strings(x):
            ''' concatenate in string and include only if it is a string (not nan), and exists '''
            return ",".join([y for y in x if isinstance(y,str) and y])
            # data = []
            # for y in x:
            #     if isinstance(y,str) and y:
            #         data.append(y)
            #     ",".join(set(data))

        def set_float_value(x):
            ''' concatenate in string and include only if it is a number, not nan, and exists '''
            return ",".join([str(y) for y in x if isinstance(y,float) and not np.isnan(y) and y])
        
        # print(self.id_xref['symbol'])

        self._search_chembl_smallmol_active()
        # from ot_tractability_pipeline_v2.buckets_ab import _process_protein_complexes
        self._process_protein_complexes()

        self.all_chembl_targets.to_csv("{}/protac_all_chembl_targets_processed.csv".format(self.store_fetched), index=False)

        # quick bug fix: only keep row when 'accession' is available
        self.all_chembl_targets = self.all_chembl_targets.loc[self.all_chembl_targets['accession'].notna()]

        self.all_chembl_targets.to_csv("{}/protac_all_chembl_targets_processed_2.csv".format(self.store_fetched), index=False)

        # self.gene_xref = self.id_xref[['accession', 'ensembl_gene_id', 'symbol']]
        # self.out_df = self.all_chembl_targets.merge(self.gene_xref, how='outer', on='accession')
        self.out_df = self.out_df.merge(self.all_chembl_targets, how='left', on='accession')
        
        # self.smallmol_active = self.out_df
        # self.smallmol_active.to_csv("{}/protac_smallmol_active.csv".format(self.store_fetched), index=False)

        # self.out_df.drop(['component_id', 'compound_name', 'ref_id', 'ref_type', 'tid', 'molregno',
        #                   'parent_molregno', 'ref_url'], axis=1, inplace=True)
        self.out_df.drop(['component_id', 'target_name', 'target_type', 'tid', #'molregno',
                          'parent_molregno'], axis=1, inplace=True)

        f = {x: 'first' for x in self.out_df.columns}
        f['compound_chembl_id'] = set_strings
        f['pchembl_value'] = set_float_value

        self.out_df = self.out_df.groupby(['ensembl_gene_id']).agg(f).reset_index(drop=True)
 
       # # quick bug fix: clear entry when 'accession' is not available
       #  self.out_df.loc[(self.out_df['accession'].isna()), 'compound_chembl_id'] = ''
       #  self.out_df.loc[(self.out_df['accession'].isna()), 'pchembl_value'] = ''

        self.out_df.rename(columns = {'compound_chembl_id':'compound_chembl_ids_PROTAC',
                                      'pchembl_value':'pchembl_values_PROTAC'}, inplace = True)

        # self.out_df.to_csv("{}/protac_out_df_checkpoint2.csv".format(self.store_fetched), index=False)

        self.out_df['Bucket_8_PROTAC'] = 0        
        # check if tuple in column 'compound_chembl_ids_PROTAC' is not empty, then set bucket to 1
        self.out_df.loc[(self.out_df['compound_chembl_ids_PROTAC'] != ''), 'Bucket_8_PROTAC'] = 1

        # Cleaning column: setting selected culumns in list format to improve visualization e.g. with Excel
        # and remove duplicates while keeping order using "list(dict.fromkeys(lst))"
        self.out_df['compound_chembl_ids_PROTAC'].fillna('', inplace=True)
        self.out_df['compound_chembl_ids_PROTAC'] = self.out_df['compound_chembl_ids_PROTAC'].apply(lambda x: list(dict.fromkeys(x.split(","))))
        # self.out_df['pchembl_values_PROTAC'].fillna('', inplace=True)
        # self.out_df['pchembl_values_PROTAC'] = self.out_df['pchembl_values_PROTAC'].apply(lambda x: list(dict.fromkeys(x.split(","))))
    
        # counting compound_chembl_ids_PROTAC: list containing empty string needs to be replaced by empty list (otherwise counted as 1)
        self.out_df['count_compound_chembl_ids_PROTAC'] = [len(x) for x in ast.literal_eval(",".join(self.out_df['compound_chembl_ids_PROTAC'].astype(str).replace(to_replace="['']",value="[]",inplace=False)))]        


        print(self.out_df.columns)


    ##############################################################################################################
    #
    # Higher level functions relating to the overall process
    #
    #
    ##############################################################################################################

    def _clinical_precedence(self, s):
        return 1 * s['Bucket_1_PROTAC'] + 0.7 * s['Bucket_2_PROTAC'] + 0.2 * s['Bucket_3_PROTAC']
    
    # def _high_conf_pred(self, s):
    #     return 0.7 * s['Bucket_4_ab'] + 0.3 * s['Bucket_5_ab']
    #
    # def _med_conf_pred(self, s):
    #     return 0.4 * s['Bucket_6_ab'] + 0.25 * s['Bucket_7_ab'] + 0.25 * s['Bucket_8_ab'] + 0.1 * s['Bucket_9_ab']

    def _summarise_buckets(self):

        print("\t- Summarising buckets...")

        self.out_df['Top_bucket_PROTAC'] = 10
        # for x in range(9, 0, -1):
        for x in range(8, 0, -1):
            self.out_df.loc[(self.out_df['Bucket_{}_PROTAC'.format(x)] == 1), 'Top_bucket_PROTAC'] = x
            self.out_df['Bucket_{}_PROTAC'.format(x)].fillna(0, inplace=True)

        self.out_df['Bucket_sum_PROTAC'] = self.out_df['Bucket_1_PROTAC'] + self.out_df['Bucket_2_PROTAC'] + \
                                           self.out_df[
                                               'Bucket_3_PROTAC'] + self.out_df['Bucket_4_PROTAC'] + self.out_df[
                                               'Bucket_5_PROTAC'] + self.out_df['Bucket_6_PROTAC'
                                           ] + self.out_df['Bucket_7_PROTAC'] + self.out_df['Bucket_8_PROTAC']# + self.out_df['Bucket_9_PROTAC']

        self.out_df.set_index('ensembl_gene_id')

    def assign_buckets(self):
        '''
        Assigns the supplied list of gene IDs into their corresponding tractability buckets.
        :return: A Pandas DataFrame containing the Ensembl gene ID and associated tractability bucket
        '''

        self._PROTAC_location_bucket()
        self._assign_buckets_1_to_3()
        self._assign_bucket_4()
        self._assign_bucket_5()
        self._assign_bucket_6()
        self._assign_bucket_7()
        self._assign_bucket_8()

        self._summarise_buckets()

        # try:

        # self.out_df.index = self.out_df['ensembl_gene_id']
        self.out_df = self.out_df.groupby('ensembl_gene_id').first()

        # Columns to keep. This includes columns from the small molecule pipeline
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
                                   'protein_name', 
                                   'Bucket_1_PROTAC', 'Bucket_2_PROTAC', 'Bucket_3_PROTAC', 
                                   'Bucket_4_PROTAC', 'Bucket_5_PROTAC', 'Bucket_6_PROTAC', 
                                   'Bucket_7_PROTAC', 'Bucket_8_PROTAC', #'Bucket_9_PROTAC', 
                                   'Bucket_sum_PROTAC', 'Top_bucket_PROTAC',
                                   # 'drug_chembl_ids_PROTAC', 'drug_names_PROTAC',
                                   'Bcell_mean', 'NKcell_mean', 'Hepatocytes_mean', 'MouseNeuorons_mean',
                                   'Max_halflife', 'Min_halflife',
                                   'Uniprot_keyword', 'Uniprot_PTM', 'Uniprot_CrossLink', 
                                   'Ub_PhosphoSitePlus', 'Ub_mUbiSiDa_2013', 'number_of_ubiquitination_sites', 
                                   'literature_count_PROTAC', 'full_id', 'title', 
                                   'count_compound_chembl_ids_PROTAC', #'compound_chembl_ids_PROTAC', #'pchembl_values_PROTAC', 
                                   'PROTAC_location_Bucket'
                                   ]]


        # Score each category, and label highest category
        self.out_df['Clinical_Precedence_PROTAC'] = self.out_df.apply(self._clinical_precedence, axis=1)
        self.out_df['Literature_Precedence_PROTAC'] = 0
        self.out_df.loc[(self.out_df['Bucket_7_PROTAC'] == 1), 'Literature_Precedence_PROTAC'] = 1
        self.out_df['Discovery_Opportunity_PROTAC'] = 0
        self.out_df.loc[(self.out_df['Bucket_4_PROTAC'] == 1) & ((self.out_df['Bucket_5_PROTAC'] == 1) | (self.out_df['Bucket_6_PROTAC'] == 1)) &
                        (self.out_df['Bucket_8_PROTAC'] == 1) & (self.out_df['PROTAC_location_Bucket'] <= 2),
                        'Discovery_Opportunity_PROTAC'] = 1

        self.out_df['Category_PROTAC'] = 'Unknown'
        self.out_df.loc[(self.out_df['Bucket_4_PROTAC'] == 1) & ((self.out_df['Bucket_5_PROTAC'] == 1) | (self.out_df['Bucket_6_PROTAC'] == 1)) &
                        (self.out_df['Bucket_8_PROTAC'] == 1) & (self.out_df['PROTAC_location_Bucket'] <= 2),
                        'Category_PROTAC'] = 'Discovery_Opportunity_PROTAC'
        self.out_df.loc[(self.out_df['Bucket_7_PROTAC'] == 1), 'Category_PROTAC'] = 'Literature_Precedence_PROTAC'
        self.out_df.loc[(self.out_df['Top_bucket_PROTAC'] <= 3), 'Category_PROTAC'] = 'Clinical_Precedence_PROTAC'
        


        print(self.out_df.columns)

        return self.out_df.astype({x: 'int64' for x in self.out_df.columns if "Bucket" in x})



    @staticmethod
    def protac2json(d):
        return {
            'Bucket_scores': {'Bucket_1_PROTAC':d.Bucket_1_PROTAC, 'Bucket_2_PROTAC':d.Bucket_2_PROTAC, 'Bucket_3_PROTAC':d.Bucket_3_PROTAC, 
                              'Bucket_4_PROTAC':d.Bucket_4_PROTAC, 'Bucket_5_PROTAC':d.Bucket_5_PROTAC, 'Bucket_6_PROTAC':d.Bucket_6_PROTAC, 
                              'Bucket_7_PROTAC':d.Bucket_7_PROTAC, 'Bucket_8_PROTAC':d.Bucket_8_PROTAC, 'PROTAC_Location_Bucket':d.PROTAC_location_Bucket}, #, 'Bucket_9_PROTAC':d.Bucket_9_PROTAC
            'Bucket_evaluation': {'Bucket_sum_PROTAC':d.Bucket_sum_PROTAC, 'Top_bucket_PROTAC':d.Top_bucket_PROTAC},
            'Category_scores': {'Clinical_Precedence_PROTAC':d.Clinical_Precedence_PROTAC, 'Literature_Precedence_PROTAC':d.Literature_Precedence_PROTAC, 'Discovery_Opportunity_PROTAC':d.Discovery_Opportunity_PROTAC}, 
            'Category_evaluation': {'Top_Category_PROTAC':d.Category_PROTAC},
            'Bucket_evidences': {'Bucket_1_2_3_PROTAC': {'drugs_in_clinic':{}, 'max_clinical_phase':{}}, #{'drug_chembl_ids_PROTAC':d.drug_chembl_ids_PROTAC}, 
                                 'Bucket_4_PROTAC': {'Max_halflife':d.Max_halflife, 'Min_halflife':d.Min_halflife}, 
                                 'Bucket_5_PROTAC': {'Uniprot_keyword':d.Uniprot_keyword, 'Uniprot_PTM':d.Uniprot_PTM, 'Uniprot_CrossLink':d.Uniprot_CrossLink}, 
                                 'Bucket_6_PROTAC': {'Ub_PhosphoSitePlus':d.Ub_PhosphoSitePlus, 'Ub_mUbiSiDa_2013':d.Ub_mUbiSiDa_2013, 'number_of_ubiquitination_sites':d.number_of_ubiquitination_sites}, 
                                 'Bucket_7_PROTAC': {'literature_count_PROTAC':d.literature_count_PROTAC, 'literature_link':d.full_id, 'title':d.title}, 
#                                 'Bucket_8_PROTAC': {'compound_chembl_ids_PROTAC':d.compound_chembl_ids_PROTAC}}
                                 'Bucket_8_PROTAC': {'number_chembl_compounds_PROTAC':d.count_compound_chembl_ids_PROTAC}}
            }
