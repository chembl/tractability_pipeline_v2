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
from sqlalchemy import create_engine # to re-establish connection
# import mygene
import numpy as np
import pandas as pd
import pkg_resources
import ast
# from sqlalchemy import create_engine

import requests

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

        # Human proteome data loaded into here
        self.human_proteome = Pipeline_setup.human_proteome
        
        # # use function from Pipeline setup
        # self.make_request = Pipeline_setup.make_request

        # All chembl data loaded into here
        self.all_chembl_targets = None

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
#            try:
#                # try reading the output from store_tetched antibody pipeline to avoid getting data again
#                uniprot_loc_data = pd.read_csv("{}/ab_uniprot_locations_processed.csv".format(self.store_fetched))
#                # self.out_uniprot_loc_data = self.out_uniprot_loc_data.merge(uniprot_loc_data, how='left', on='accession')
#
#            except IOError:
#                print("--- \nWarning in PROTAC workflow: Uniprot location file 'ab_uniprot_locations_processed.csv' is not available \
#                      from antibody workflow. (see section 'Uniprot (loc)' in buckets_ab.py) \n-> trying to reasses from source \nThis will take some time ...\n---")
#                uniprot_loc_data = None
#            
#            if uniprot_loc_data is None:
#                #=============================================================================
#                #    The following code is only executed if no previously calculated ab_uniprot_locations_processed.csv is found     
#                #=============================================================================
#                
#                from ot_tractability_pipeline_v2.buckets_ab import Antibody_buckets
#                # make_request = Antibody_buckets.make_request
#                # post_request_uniprot = Antibody_buckets.post_request_uniprot
#                # split_loc = Antibody_buckets.split_loc
#                # _check_evidence = Antibody_buckets._check_evidence
#                            
#                # Return all reviewed and Human targets
#                url = "uniprot/?format=tab&query=*&fil=reviewed%3ayes+AND+organism%3a%22Homo+sapiens+(Human)+%5b9606%5d%22&columns=id,comment(SUBCELLULAR+LOCATION),comment(DOMAIN),feature(DOMAIN+EXTENT),feature(INTRAMEMBRANE),feature(TOPOLOGICAL+DOMAIN),feature(TRANSMEMBRANE),feature(SIGNAL)"
#                data = ['P42336', 'P60484']
#                location = Antibody_buckets.post_request_uniprot(url, data)
#                location = [x.split('\t') for x in location.split('\n')]
#                uniprot_loc_data = pd.DataFrame(location[1:], columns=location[0])
#                uniprot_loc_data['uniprot_loc_test'] = uniprot_loc_data['Subcellular location [CC]']
#                uniprot_loc_data['Subcellular location [CC]'] = uniprot_loc_data['Subcellular location [CC]'].apply(Antibody_buckets.split_loc)
#        
#                uniprot_loc_data.rename(columns={'Entry': 'accession'}, inplace=True)

        if 'Uniprot_high_conf_loc' not in self.out_df.columns:
            # =============================================================================
            #     The following code is only executed if no previously generated Uniprot_high_conf_loc column is found     
            # =============================================================================
            print("\t--- Processed UniProt location data (columns 'Uniprot_high_conf_loc' and 'Uniprot_med_conf_loc') is not available. \
                  Processing UniProt location data ... ---")
            from ot_tractability_pipeline_v2.buckets_ab import Antibody_buckets
            # UniProt location data is already fetched in Pipeline setup
            # Copy original data column for eventual later comparison
            self.out_df['original_uniprot_location_data'] = self.out_df['Subcellular location [CC]']
            # Process the delimited string returned from uniprot webservice call
            self.out_df['Subcellular location [CC]'] = self.out_df['Subcellular location [CC]'].apply(Antibody_buckets.split_loc)
            # Assign UniProt high and medium confidence locations 

            def get_uniprot_high_conf(s):            
                all_uniprot_high_conf = [(a[1], a[0]) for a in s['Subcellular location [CC]'] if Antibody_buckets._check_evidence(a[0])]    
                return all_uniprot_high_conf

            def get_uniprot_med_conf(s):            
                all_uniprot_med_conf = [(a[1], a[0]) for a in s['Subcellular location [CC]'] if not Antibody_buckets._check_evidence(a[0])]
                return all_uniprot_med_conf
    
            self.out_df['Uniprot_high_conf_loc'] = self.out_df.apply(get_uniprot_high_conf, axis=1)
            self.out_df['Uniprot_med_conf_loc'] = self.out_df.apply(get_uniprot_med_conf, axis=1)

#                if self.store_fetched: 
#                    uniprot_loc_data.to_csv("{}/protac_uniprot_locations_processed.csv".format(self.store_fetched))

        else:
            print("\t--- Already available processed UniProt location data is used. ---")

        # finally:
        #print(uniprot_loc_data.columns)
        #self.out_df = self.out_df.merge(uniprot_loc_data, how='left', on='accession')
        print(self.out_df.columns)

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
   
        if go_loc_data is None:
        # =============================================================================
        #     The following code is only executed if no previously calculated ab_GO_locations_processed.csv is found     
        # =============================================================================
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
#                    accepted_high_conf_loc = []
#                    accepted_med_conf_loc = []
        
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
        
                return high_conf_loc, med_conf_loc
        
            '''
            GO CC
            '''
            self.out_df['GO_high_conf_loc'], self.out_df['GO_med_conf_loc'] = zip(*self.out_df.apply(get_go_loc, axis=1))
            
            # save processed GO location info to file
            if self.store_fetched: 
                self.out_df.loc[:, self.out_df.columns.isin(['accession', 'GO_high_conf_loc', 'GO_med_conf_loc'])].to_csv(
                    "{}/protac_GO_locations_processed.csv".format(self.store_fetched))
    
            print(self.out_df.columns)
            

        else:
            # print("\t--- GO location data used from Antibody workflow ---")
            #print(go_loc_data.columns)
            print(self.out_df.columns)


    # ChEMBL currently not used

    # if database_url is None:
    #     database_url = os.getenv('CHEMBL_DB')
    #
    #
    # # Create ChEMBL DB connection
    # self.engine = create_engine(database_url)        


    ##############################################################################################################
    #
    # Functions relating to PROTAC_location_bucket
    # 
    #
    ##############################################################################################################


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
        # elif len(bad_locations) > 0:
        #     return 6
        # If medium conf locations are known, but not in self.good_locations or self.grey_locations, they are assumed to be bad
        # else:
        #     return 6
        elif len(locations) > 0:
            return 6
 

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
            # If information on GO CC is provided by uniprot return as medium confidence good/grey location
            if row['Gene ontology (cellular component)'] != '' and any(x in row['Gene ontology (cellular component)'] for x in self.good_locations):
                return 2
            elif row['Gene ontology (cellular component)'] != '' and any(x in row['Gene ontology (cellular component)'] for x in self.grey_locations):
                return 4
            else: 
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
        # elif len(bad_locations) > 0:
        #     return 7
        
        # Taking into account medium conf assessment (not overriding score 2 and 4; only 6); If high conf locations are known, but not in self.good_locations or self.grey_locations, they are assumed to be bad
        elif (row['PROTAC_location_Bucket'] == 6) and (len(locations) > 0):
            return 7
#        # If high conf locations are known, but not in self.good_locations or self.grey_locations, they are assumed to be bad
#        elif len(locations) > 0:
#            return 7
        else:
            #print('locations',locations, 'accepted', accepted_locations, 'grey', grey_locations)
            return row['PROTAC_location_Bucket']


    def _PROTAC_location_bucket(self):

        ''''
        For PROTACs, only intracellular targets are suitable. Therefore, we will assign a score based on location to
        allow the targets to be filtered to those in the cytosol, nucleus or membrane with accessible portion.

        1 - High confidence good location
        2 - Med confidence good location
        3 - High confidence grey location
        4 - Med condfidence grey location
        5 - Unknown location
        6 - Med confidence bad location
        7 - High confidence bad location
        '''

        print("\t- Assessing PROTAC location bucket...")

        self.good_locations = ['cytoplasm', 'cytosol', 'nucleus', 'nucleolus', 'nucleoplasm', 'nucleus matrix']
        self.grey_locations = ["membrane", "nuclear pore complex", 
                               "cell membrane", "plasma membrane", "cytoplasmic membrane", "plasmalemma", 
                               "apical cell membrane", "apical plasma membrane", 
                               "apicolateral cell membrane", "apicolateral plasma membrane", 
                               "basal cell membrane", "basal plasma membrane", 
                               "basolateral cell membrane", "basolateral plasma membrane", 
                               "lateral cell membrane", "lateral plasma membrane", 
                               "synaptic cell membrane", "presynaptic cell membrane", "postsynaptic cell membrane", 
                               "synaptic vesicle membrane", "pseudopodium membrane", "bud membrane", "caveola", 
                               "filopodium membrane", "flagellum membrane", "cilium membrane", "ciliary membrane",
                               "endoplasmic reticulum membrane", "rough endoplasmic reticulum membrane", "smooth endoplasmic reticulum membrane", 
                               "Golgi membrane", "golgi membrane", "golgi apparatus membrane", "golgi stack membrane", "trans-golgi network membrane",
                               "peroxisome membrane", "endosome membrane", "phagosome membrane", "esterosome membrane", "lysosome membrane", 
                               "mitochondrion outer membrane", "nucleus envelope", "nucleus membrane", 
                               "attachment organelle membrane", "cellular thylakoid membrane", "plastid outer membrane", "vacuole membrane"]
        # self.bad_locations = ['secreted']

        self.out_df['PROTAC_location_Bucket'] = 0
        self.out_df['PROTAC_location_Bucket'] = self.out_df.apply(self._med_conf_locations, axis=1)
        self.out_df['PROTAC_location_Bucket'] = self.out_df.apply(self._high_conf_locations, axis=1)

        print(self.out_df.columns)
        if self.store_fetched: 
            self.out_df.to_csv("{}/Checkpoint_PROTAC.csv".format(self.store_fetched))



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
        #self.out_df.loc[(self.out_df['ensembl_gene_id'] == 'ENSG00000169083'), 'Bucket_3_PROTAC'] = 1
        self.out_df.loc[(self.out_df['ensembl_gene_id'] == 'ENSG00000169083'), 'Bucket_2_PROTAC'] = 1
        # ER PROTAC ARV-471 is currently in clinical trials for advanced or metastatic ER+/Her2− breast cancer. 
        #self.out_df.loc[(self.out_df['ensembl_gene_id'] == 'ENSG00000091831'), 'Bucket_3_PROTAC'] = 1
        self.out_df.loc[(self.out_df['ensembl_gene_id'] == 'ENSG00000091831'), 'Bucket_2_PROTAC'] = 1
        # IRAK4 PROTAC KT-474 is currently in clinical trials for the treatment of interleukin-1 receptor (IL-1R)/toll-like receptor (TLR)-driven immune-inflammatory diseases.
        self.out_df.loc[(self.out_df['ensembl_gene_id'] == 'ENSG00000198001'), 'Bucket_3_PROTAC'] = 1
        # Bruton’s tyrosine kinase (BTK) PROTAC NX-2127 is currently in clinical trials for the treatment of relapsed or refractory B-cell malignancies.
        self.out_df.loc[(self.out_df['ensembl_gene_id'] == 'ENSG00000010671'), 'Bucket_3_PROTAC'] = 1

        # Cleaning column: setting selected culumns in list format to improve visualization e.g. with Excel
        # and remove duplicates while keeping order using "list(dict.fromkeys(lst))"
        # self.out_df['drug_chembl_ids_PROTAC'].fillna('', inplace=True)
        # self.out_df['drug_chembl_ids_PROTAC'] = self.out_df['drug_chembl_ids_PROTAC'].apply(lambda x: list(dict.fromkeys(x.split(","))))

        print(self.out_df.columns)


    ##############################################################################################################
    #
    # Functions relating to bucket 4(+)
    # Manually curated PROTAC targets +
    # Additional info: Taregts mentioned in PROTAC literature, automatically detected
    #
    ##############################################################################################################

    @staticmethod
    def _search_papers():
        
        # keyword search only in Abstract (keywords in "" need additional 'ABSTRACT:'): 
        # ABSTRACT:((ABSTRACT:"proteolysis targeting chimera" OR ABSTRACT:"proteolysis targeting chimeric") OR ((PROTAC OR SNIPER) AND (degradation OR degrade OR proteolysis)) OR (degrader AND (proteasome OR ubiquitin)))
        # use resultType=core& before query to get full abstract text in return
        # ABSTRACT:((ABSTRACT:"proteolysis targeting chimera" OR ABSTRACT:"proteolysis targeting chimeric") OR (PROTAC AND (degradation OR degrade OR ubiquitin OR proteolysis)))
        #url = "https://www.ebi.ac.uk/europepmc/webservices/rest/search?resultType=core&query=ABSTRACT%3A%28%28ABSTRACT%3A%22proteolysis%20targeting%20chimera%22%20OR%20ABSTRACT%3A%22proteolysis%20targeting%20chimeric%22%29%20OR%20%28PROTAC%20AND%20%28degradation%20OR%20degrade%20OR%20ubiquitin%20OR%20proteolysis%29%29%29&resultType=lite&cursorMark=*&pageSize=1000&format=json"
    
        # urllib not working returning 404 page error -> replaced by requests
        #url = urllib2.urlopen("https://www.ebi.ac.uk/europepmc/webservices/rest/search?resultType=core&query=ABSTRACT%3A%28%28ABSTRACT%3A%22proteolysis%20targeting%20chimera%22%20OR%20ABSTRACT%3A%22proteolysis%20targeting%20chimeric%22%29%20OR%20%28%28PROTAC%20OR%20SNIPER%29%20AND%20%28degradation%20OR%20degrade%20OR%20proteolysis%29%29%20OR%20%28degrader%20AND%20%28proteasome%20OR%20ubiquitin%29%29%29&resultType=lite&cursorMark=*&pageSize=1000&format=json")
        #data = url.read()
        #try: data = json.loads(data.decode())
        #except UnicodeDecodeError: data = json.loads(data)

        url = "https://www.ebi.ac.uk/europepmc/webservices/rest/search?resultType=core&query=ABSTRACT%3A%28%28ABSTRACT%3A%22proteolysis%20targeting%20chimera%22%20OR%20ABSTRACT%3A%22proteolysis%20targeting%20chimeric%22%29%20OR%20%28%28PROTAC%20OR%20SNIPER%29%20AND%20%28degradation%20OR%20degrade%20OR%20proteolysis%29%29%20OR%20%28degrader%20AND%20%28proteasome%20OR%20ubiquitin%29%29%29&resultType=lite&cursorMark=*&pageSize=1000&format=json"
        response = requests.get(url)
        data = response.json()
        
        #df = pd.read_json(json.dumps(data['resultList']['result']), orient='records')
        df = pd.json_normalize(data['resultList']['result'])
    
        # remove entries with missing 'authorString'
        df = df.loc[~df['authorString'].isna()]
        
        #"['journalTitle', 'issue', 'pubType'] not in index"
        return df[['authorString', 'id', 'pmcid',
                   'pmid', 'pubYear', 'source', 'title', 'abstractText']]
    

    @staticmethod
    def _search_ID(row):
        return "articleIds={}%3A{}".format(row['source'], row['id'])

    @staticmethod
    def _full_ID(row):
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

    
# =============================================================================
#     @staticmethod
#     def _get_human_proteome():
#         '''
#         Getting the human proteome from UniProt; processing entry, gene, and protein names
#         '''
#         full_url = 'https://www.uniprot.org/uniprot/?query=proteome:UP000005640&format=tab&columns=id,entry%20name,protein%20names,genes'
#         from ot_tractability_pipeline_v2.buckets_ab import Antibody_buckets
#         Uniprot_human_proteome = Antibody_buckets.make_request(full_url, data=None)
#         Uniprot_human_proteome = [x.split('\t') for x in Uniprot_human_proteome.split('\n')]
#         human_proteome = pd.DataFrame(Uniprot_human_proteome[1:], columns=Uniprot_human_proteome[0])
#         human_proteome.rename(columns={'Entry': 'accession'}, inplace=True)
#         # only keep row when 'Entry name' is available (discard NAN row)
#         human_proteome = human_proteome.loc[human_proteome['Entry name'].notna()]
#         # create 'symbol' column
#         human_proteome[['symbol','Human']] = human_proteome['Entry name'].str.split("_",expand=True)
#         # create 'gene_name' column (using first entry in 'Gene names')
#         human_proteome['gene_name'] = human_proteome['Gene names'].str.split(" ",expand=True)[0]
#         # create 'protein_name' column (using primary entry, before names in parentheses in 'Protein names', escaping with \\ is reqired)
#         human_proteome['protein_name'] = human_proteome['Protein names'].str.split(" \\(",expand=True)[0]
#         # save lower case protein name for future case insensitive mapping
#         human_proteome['protein_name_lower'] = human_proteome['protein_name'].str.lower()
#         # as all protein isoforms have different UniProt IDs, only the first occurence of gene_name is kept 
#         # (which should be the primary UniProtID) count: 20487
#         human_proteome.drop_duplicates(subset="gene_name", keep='first', inplace=True)
#         
#         return human_proteome
# 
# =============================================================================

    def _process_IDs(self):
        '''
        1) Processing the EuropePMC literature ids, the provided annotations and the mappings to UniProt accessions.
        2) Getting the human proteome from Uniprot in order to recover more mappings by 'symbol', 'gene_name' and 'protein_name'.
        '''
        # 1) ---
        #grouped_tags = self.tags.groupby('name').first()
        #tagged_annotations = self.annotations.merge(grouped_tags, how='left', left_on='exact', right_on='name')
        tagged_annotations = self.annotations.merge(self.tags, how='left', left_index=True, right_index=True)

        # keep only tagged_annotations from Abstract of publications
        tagged_annotations = tagged_annotations.loc[tagged_annotations['section'].astype(str).str.contains('Abstract')]
        
        tagged_annotations['EuropePMC_accession'] = tagged_annotations.apply(self._extract_uniprot, axis=1)
        
        tagged_annotations['short_id'] = tagged_annotations.apply(self._extract_id, axis=1)
        
        tagged_annotations.rename(columns={'id':'link_id'}, inplace=True)
        
        # join tagged_annotations with self.papers_df
        joined = tagged_annotations.merge(self.papers_df, left_on='short_id', right_on='id', how='inner')
        # convert name and protein_name to lowercase for merging 
        #(to avoid missing through different uppercase usage in text)
        joined['name_lower'] = joined['name'].str.lower()
        
        # 2) ---
        #human_proteome = self._get_human_proteome()
        human_proteome = self.human_proteome
        # Join on human_proteome data to get human UniProt accession IDs:
        # 1. join with human_proteome on 'symbol' (higher chances of getting the correct/dominant isoform)
        # 2. join with human_proteome on 'gene_name' (occurs several times, for all isoforms)
        # 3. join both new datasets
        tagged_targets_on_symbol = joined.merge(human_proteome, left_on='name', right_on='symbol', how='left')
        tagged_targets_on_gene_name = joined.merge(human_proteome, left_on='name', right_on='gene_name', how='left')
        tagged_targets_df = tagged_targets_on_symbol.merge(tagged_targets_on_gene_name[['name','accession','gene_name']], on='name', how='inner')
        # 4. join with human_proteome on 'protein_name_lower' (to get full name annotations) = third new dataset 'tagged_targets_on_protein_name'
        # 5. join tagged_targets_df with third new dataset 'tagged_targets_on_protein_name' using lowercase names
        tagged_targets_on_protein_name = joined.merge(human_proteome, left_on='name_lower', right_on='protein_name_lower', how='left')
        tagged_targets_df = tagged_targets_df.merge(tagged_targets_on_protein_name[['name','accession','protein_name']], on='name', how='left')
        # 6. remove duplicated rows (keep='first' is default), as lists are contained in df:
        # convert the df to str type astype(str), drop duplicates and then select the rows from original df, thus output df still contains lists
        tagged_targets_df = tagged_targets_df.loc[tagged_targets_df.astype(str).drop_duplicates().index]        
        # 7. rename last joined accession column to "accession_z"
        tagged_targets_df.rename(columns={"accession": "accession_z"}, inplace=True)        
        # 8. add additional mapping 'accession_z2' column based on protein name from publication matching end of UniProts full protein name
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

        # safe human_proteome into global variable for later use in other buckets
        #self.human_proteome = human_proteome[['accession','symbol','gene_name','protein_name']]

        if self.store_fetched: 
            tagged_targets_df.to_csv("{}/protac_pmc_tagged_targets_raw.csv".format(self.store_fetched), encoding='utf-8')
        
        tagged_targets_df.rename(columns = {'id':'pub_id'}, inplace = True)
#        return joined[['accession', 'prefix', 'exact', 'postfix', 'section', 'full_id', 'journalTitle', 'title']]
        return tagged_targets_df[['accession', 'name', 'exact', 'pub_id', 'short_id', 'full_id', 'section', 'title', 
                                  'abstractText', 'NER_PROTAC_label']] #, 'sentences', 'sentence_count'

     
    def _detect_targets_in_literature(self, tagged_targets_df):
        '''
        Reformat a copy of tagged_targets to provide the exact tagged term and to be mergable by paper link and UniProt accession. 
        '''
        # get paper links and tagged proteins
        tagged_targets_links = tagged_targets_df.drop_duplicates(['exact','full_id','section'],ignore_index=True)[['exact','accession','full_id','section']]
        # keep only annotations from Abstracts, then drop column 'section'
        tagged_targets_links = tagged_targets_links[tagged_targets_links['section'].str.contains("Abstract")].drop(columns=['section'])
        # group by paper and set strings for protein search terms
        def set_strings(x):
            ''' concatenate in string and include only if it is a string (not nan), and exists '''
            return ", ".join([y for y in x if isinstance(y,str) and y])
        
        f0 = {x: set_strings for x in tagged_targets_links.columns if x not in ['accession','full_id']}
        
        tagged_targets_links = tagged_targets_links.groupby(['full_id','accession'], as_index=False).agg(f0).reset_index(drop=True)
        tagged_targets_links.rename(columns = {'exact':'exact_terms'}, inplace = True)
        
        # Merge tags with paper data
        tagged_targets_abstract = tagged_targets_links.merge(self.papers_df, how='left', on='full_id') #.drop(columns=['abstractText','sentences','sentence_count'])
        
        return tagged_targets_abstract
    

    @staticmethod
    def _get_human_ubi_ligases():
        '''
        Getting the human proteins taking part in ubiquitin ligase complex from UniProt using curated function and interaction fields (avoiding degradation targets from Post-translational modification field)
        OLD query: ("ubiquitin ligase complex" OR name:"ubiquitin-protein ligase") AND organism:"Homo sapiens (Human) [9606]" AND proteome:up000005640
        query: (annotation:((type:function "ubiquitin-protein ligase complex") OR (type:function "ubiquitin ligase complex") OR (type:function "e3 protein ligase") OR (type:subunit "ubiquitin-protein ligase") OR (type:subunit "ubiquitin ligase complex")) OR name:"ubiquitin-protein ligase") AND organism:"Homo sapiens (Human) [9606]" AND proteome:up000005640
        '''
        full_url = 'https://www.uniprot.org/uniprot/?query=%28annotation%3A%28%28type%3Afunction+%22ubiquitin-protein+ligase+complex%22%29+OR+%28type%3Afunction+%22ubiquitin+ligase+complex%22%29+OR+%28type%3Afunction+%22e3+protein+ligase%22%29+OR+%28type%3Asubunit+%22ubiquitin-protein+ligase%22%29+OR+%28type%3Asubunit+%22ubiquitin+ligase+complex%22%29%29+OR+name%3A%22ubiquitin-protein+ligase%22%29+AND+organism%3A%22Homo+sapiens+%28Human%29+%5B9606%5D%22+AND+proteome%3Aup000005640&format=tab&columns=id,entry+name,reviewed,protein+names,genes'
        from ot_tractability_pipeline_v2.buckets_ab import Antibody_buckets
        Uniprot_human_ubi_ligases = Antibody_buckets.make_request(full_url, data=None)
        Uniprot_human_ubi_ligases = [x.split('\t') for x in Uniprot_human_ubi_ligases.split('\n')]
        human_ubi_ligases = pd.DataFrame(Uniprot_human_ubi_ligases[1:], columns=Uniprot_human_ubi_ligases[0])
        human_ubi_ligases.rename(columns={'Entry': 'accession'}, inplace=True)
        # only keep row when 'Entry name' is available (discard NAN row)
        human_ubi_ligases = human_ubi_ligases.loc[human_ubi_ligases['Entry name'].notna()]
#        # create 'gene_name' column (using first entry in 'Gene names')
#        human_ubi_ligases['gene_name'] = human_ubi_ligases['Gene names'].str.split(" ",expand=True)[0]
#        # as all protein isoforms have different UniProt IDs, only the first occurence of gene_name is kept 
#        # (which should be the primary UniProtID) count: 20487
#        human_ubi_ligases.drop_duplicates(subset="gene_name", keep='first', inplace=True)
        
        return human_ubi_ligases


    def _assign_bucket_4(self):
        '''
        Manually curated PROTAC literature
        '''
        print("\t- Assessing PROTAC literature bucket 4...")
        
        protac_targets = pd.read_csv(os.path.join(DATA_PATH, 'PROTAC_targets.csv'))
        self.out_df = self.out_df.merge(protac_targets[['accession','name_used','protac_target']], how='left', on='accession')
        self.out_df = self.out_df.dropna(subset=['accession'])
        
        self.out_df['Bucket_4_PROTAC'] = 0
        self.out_df.loc[(~self.out_df['protac_target'].isna()), 'Bucket_4_PROTAC'] = 1
        
        '''
        Degradable kinome from Fischer et al 2020
        '''
        deg_kinome = pd.read_excel(os.path.join(DATA_PATH, 'degradable_kinome_Fischer2020_S4_mmc5.xlsx'), sheet_name='All - freq isoform collapse ')
        #self.out_df = self.out_df.merge(deg_kinome, how='left', left_on='accession', right_on='Accession')
        self.out_df = self.out_df.merge(deg_kinome, how='left', left_on='symbol', right_on='ID')
        self.out_df.loc[(~self.out_df['Freq.down'].isna()), 'Bucket_4_PROTAC'] = 1
        
        
        '''
        Mentioned in PROTAC literature (automated detection)
        '''
        print("\t- Assessing PROTAC literature (automated detection)...")

        self.papers_df = self._search_papers()
        
        self.papers_df['search_id'] = self.papers_df.apply(self._search_ID, axis=1)
        self.papers_df['full_id'] = self.papers_df.apply(self._full_ID, axis=1)

        # =============================================================================        
        # Use pretrained SpaCy model for Named Entity Recognition trained for detection of 'PROTAC_TARGET', 'PROTAC_NAME', and 'E3_LIGASE'.

        import spacy
#        model_dir = os.path.join(DATA_PATH, "SpaCy_NER_PROTAC_model")
#        print("Loading SpaCy model for Named Entity Recognition of 'PROTAC_TARGET', 'PROTAC_NAME', and 'E3_LIGASE' from", model_dir)
#        nlp2 = spacy.load(model_dir)
        print("Loading SpaCy model for Named Entity Recognition of 'PROTAC_TARGET', 'PROTAC_NAME', and 'E3_LIGASE'.")
        nlp2 = spacy.load(r'en_NER_PROTAC')
        
        def predict_abstracts(textcol):
            docx = nlp2(textcol)
            pred = [(ent.label_, ent.text) for ent in docx.ents]
            return pred
                
        self.papers_df['NER_PROTAC_label'] = self.papers_df['abstractText'].map(predict_abstracts)

        # =============================================================================

        if self.store_fetched: 
            self.papers_df.to_csv("{}/protac_pmc_papers.csv".format(self.store_fetched), encoding='utf-8')
                
        self._get_tagged_targets()

        tagged_targets_df = self._process_IDs()
        
        tagged_targets_df = tagged_targets_df.drop_duplicates(['accession','full_id'], ignore_index=True)
        #tagged_targets_df = tagged_targets_df.loc[tagged_targets_df.astype(str).drop_duplicates().index]

        # remove duplicated rows (keep='first' is default), as lists are contained in df:
        # convert the df to str type astype(str), drop duplicates and then select the rows from original df, thus output df still contains lists
        tagged_targets_df = tagged_targets_df.loc[tagged_targets_df.astype(str).drop_duplicates(ignore_index=True).index]

        def set_strings(x):
            ''' concatenate in string and include only if it is a string (not nan), and exists '''
            return ",".join([y for y in x if isinstance(y,str) and y])
        def set_title_strings(x):
            ''' concatenate in string and include only if it is a string (not nan), and exists '''
            return " | ".join([y for y in x if isinstance(y,str) and y])
        def set_list(x):
            return list(x)
        def set_from_list(x):
            ''' concatenate in string and include only if it is a list and not empty '''
            return " | ".join([str(y) for y in x if isinstance(y,list) and len(y)>0])
            #apply(list) #if isinstance(a,list) and len(a) > 0
        
        f = {x: 'first' for x in tagged_targets_df.columns} # if x != 'accession'
        f['pub_id'] = set_strings
        f['short_id'] = set_strings
        f['full_id'] = set_strings
        f['title'] = set_title_strings
        #f['exact_terms'] = set_title_strings
        f['abstractText'] = set_title_strings
        f['NER_PROTAC_label'] = set_from_list
        
        tagged_targets_df_grouped = tagged_targets_df.groupby(['accession']).agg(f).reset_index(drop=True)
        
        # removing duplicates        
        def clean_column(df, col):
            df[col] = df[col].apply(lambda x: list(dict.fromkeys(x.split(","))))
        def clean_title_column(df, col):
            df[col] = df[col].apply(lambda x: list(dict.fromkeys(x.split(" | "))))
        
        clean_column(tagged_targets_df_grouped, 'pub_id')
        clean_column(tagged_targets_df_grouped, 'short_id')
        clean_column(tagged_targets_df_grouped, 'full_id')
        clean_title_column(tagged_targets_df_grouped, 'title')
        clean_title_column(tagged_targets_df_grouped, 'abstractText')
        #clean_title_column(tagged_targets_df_grouped, 'sentences')
        
        tagged_targets_df_grouped = tagged_targets_df_grouped.loc[tagged_targets_df_grouped.astype(str).drop_duplicates(ignore_index=True).index]
        if self.store_fetched: 
            tagged_targets_df.to_csv("{}/protac_pmc_tagged_targets_processed.csv".format(self.store_fetched), encoding='utf-8')

        self.out_df = self.out_df.merge(tagged_targets_df_grouped, how='left', on='accession')
        
        self.out_df['mentioned_in_PROTAC_literature'] = 0
        self.out_df.loc[(~self.out_df['full_id'].isna()), 'mentioned_in_PROTAC_literature'] = 1
        
        # count literature per target and append as column
        self.out_df['literature_count_PROTAC'] = [len(x) for x in ast.literal_eval(",".join(self.out_df['full_id'].fillna("['']",inplace=False).astype(str).replace(to_replace="['']",value="''",inplace=False)))]

        ## adding protein_name from human_proteome to self.out_df
        #self.out_df = self.human_proteome[['accession','protein_name']].merge(self.out_df, how='right', on='accession')        
        # removing duplicated rows        
        self.out_df = self.out_df.drop_duplicates(['ensembl_gene_id'], ignore_index=True)
        #self.out_df = self.out_df.loc[self.out_df.astype(str).drop_duplicates(ignore_index=True).index]
        
        # save additional evaluation results
        if self.store_fetched: 
            literature_evidence = self.out_df.loc[self.out_df['mentioned_in_PROTAC_literature'] == 1]
            literature_evidence = literature_evidence[['accession','symbol','protein_name','literature_count_PROTAC','title','full_id',
                                 'abstractText', 'NER_PROTAC_label']].sort_values(
                    by=['literature_count_PROTAC'],ascending=[0]) # ,'sentences','sentence_count','tag_in_sentence_count'
            literature_evidence = literature_evidence.drop_duplicates(['accession'], ignore_index=True)
            literature_evidence.to_csv("{}/protac_literature_evidence.csv".format(self.store_fetched), encoding='utf-8')

            human_ubi_ligases = self._get_human_ubi_ligases()
            literature_evidence_ubi_ligases = literature_evidence.merge(human_ubi_ligases['accession'], how='inner')
            literature_evidence_ubi_ligases.to_csv("{}/protac_literature_evidence_ubi_ligase_complex.csv".format(self.store_fetched), encoding='utf-8')
            #literature_evidence_ubi_ligases = literature_evidence[literature_evidence['accession'].isin(human_ubi_ligases['accession'])]
            #literature_evidence_ubi_ligases.to_csv("{}/protac_literature_evidence_ubi_ligases2.csv".format(self.store_fetched), encoding='utf-8')
            literature_evidence_filtered = literature_evidence[~literature_evidence['accession'].isin(human_ubi_ligases['accession'])]
            literature_evidence_filtered.to_csv("{}/protac_literature_evidence_filtered.csv".format(self.store_fetched), encoding='utf-8')

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
        #self.out_df = ub_df.merge(self.out_df, on='symbol', how='right')
        self.out_df = self.out_df.merge(ub_df, on='symbol', how='left')

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
    # Protein Turnover
    #
    ##############################################################################################################

    def _assign_bucket_7(self):
        '''
        Protein Turnover
        '''

        print("\t- Assessing protein turnover bucket 7...")

        self.out_df['Bucket_7_PROTAC'] = 0
        # self.out_df['Bucket_x_PROTAC'] = 0

        df = pd.read_csv(os.path.join(DATA_PATH, 'protein_half_life_hq.csv'))

        df = df.merge(self.out_df, right_on='symbol', left_on='gene_name', how='right')
        df = df.groupby('ensembl_gene_id', as_index=False).max(numeric_only=True)
        df['Max_halflife'].fillna(-1, inplace=True)
        df['Min_halflife'].fillna(-1, inplace=True)

        self.out_df = df.merge(self.out_df, how='right', on='ensembl_gene_id', suffixes=['_drop', ''])
        
        # Assign bucket if Half-life data is available
        self.out_df.loc[(self.out_df['Min_halflife'] > 0), 'Bucket_7_PROTAC'] = 1
        # self.out_df.loc[(self.out_df['Max_halflife'] >= 24), 'Bucket_7_PROTAC'] = 1
        # self.out_df.loc[(self.out_df['Max_halflife'] > 10) & (self.out_df['Max_halflife'] < 24), 'Bucket_x_PROTAC'] = 1

        print(self.out_df.columns)



    ##############################################################################################################
    #
    # Functions relating to bucket 8
    # Existing active molecule
    # (before: Small Molecule Tractable)
    #
    ##############################################################################################################

    # def _assign_bucket_8(self):
    #     '''
    #     Small molecule tractable
    #     '''
    #     self.out_df['Bucket_9_PROTAC'] = 0
    #     self.out_df.loc[(self.out_df['Top_bucket_sm'] < 9), 'Bucket_9_PROTAC'] = 1


    def _search_chembl_smallmol_active(self):
        '''
        Search for all targets in ChEMBL, get all small molecules and their targets with a pChEMBL >=5
        (see queries_protac.py)
        '''
        # print("\t\t- Querying ChEMBL...")
        
        # URL to local ChEMBL database
        # If database url not supplied, get from environemnt variable
        database_url = None
        if database_url is None:
            database_url = os.getenv('CHEMBL_DB')
        # ChEMBL DB connection - reestablish connection to avoid timeout
        self.engine = create_engine(database_url)

        # small_mol_info = pd.read_sql_query(chembl_small_mol, self.engine)
        self.all_chembl_targets = pd.read_sql_query(chembl_small_mol_active_targets, self.engine)
        # self.all_chembl_targets = self.all_chembl_targets.merge(small_mol_info, on='parent_molregno')

        if self.store_fetched: 
            self.all_chembl_targets.to_csv("{}/protac_all_chembl_active_targets.csv".format(self.store_fetched))


    def _process_protein_complexes(self):
        '''
        For protein complexes, see if we know the binding subunit, and only keep these
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
        '''
        print("\t- Assessing existing active molecule bucket 8...")

        def set_strings(x):
            ''' concatenate in string and include only if it is a string (not nan), and exists '''
            return ",".join([y for y in x if isinstance(y,str) and y])

        def set_float_value(x):
            ''' concatenate in string and include only if it is a number, not nan, and exists '''
            return ",".join([str(y) for y in x if isinstance(y,float) and not np.isnan(y) and y])
        
        self._search_chembl_smallmol_active()
        # from ot_tractability_pipeline_v2.buckets_ab import _process_protein_complexes
        self._process_protein_complexes()

        self.all_chembl_targets.to_csv("{}/protac_all_chembl_targets_processed.csv".format(self.store_fetched), index=False)

        # quick bug fix: only keep row when 'accession' is available
        self.all_chembl_targets = self.all_chembl_targets.loc[self.all_chembl_targets['accession'].notna()]

        # self.gene_xref = self.id_xref[['accession', 'ensembl_gene_id', 'symbol']]
        # self.out_df = self.all_chembl_targets.merge(self.gene_xref, how='outer', on='accession')
        self.out_df = self.out_df.merge(self.all_chembl_targets, how='left', on='accession')
        
        # self.out_df.drop(['component_id', 'compound_name', 'ref_id', 'ref_type', 'tid', 'molregno',
        #                   'parent_molregno', 'ref_url'], axis=1, inplace=True)
        self.out_df.drop(['component_id', 'target_name', 'target_type', 'tid', #'molregno',
                          'parent_molregno'], axis=1, inplace=True)

        f = {x: 'first' for x in self.out_df.columns}
        f['compound_chembl_id'] = set_strings
        f['pchembl_value'] = set_float_value

        self.out_df = self.out_df.groupby(['ensembl_gene_id']).agg(f).reset_index(drop=True)
 
        # # quick bug fix: clear entry when 'accession' is not available (not needed any more)
        # self.out_df.loc[(self.out_df['accession'].isna()), 'compound_chembl_id'] = ''
        # self.out_df.loc[(self.out_df['accession'].isna()), 'pchembl_value'] = ''

        self.out_df.rename(columns = {'compound_chembl_id':'compound_chembl_ids_PROTAC',
                                      'pchembl_value':'pchembl_values_PROTAC'}, inplace = True)

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

        # self.out_df.index = self.out_df['ensembl_gene_id']
        self.out_df = self.out_df.groupby('ensembl_gene_id').first()

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
                                   'Bucket_1_PROTAC', 'Bucket_2_PROTAC', 'Bucket_3_PROTAC', 
                                   'Bucket_4_PROTAC', 'Bucket_5_PROTAC', 'Bucket_6_PROTAC', 
                                   'Bucket_7_PROTAC', 'Bucket_8_PROTAC', #'Bucket_9_PROTAC', 
                                   'Bucket_sum_PROTAC', 'Top_bucket_PROTAC',
                                   'PROTAC_location_Bucket',
                                   # 'drug_chembl_ids_PROTAC', 'drug_names_PROTAC',
                                   'mentioned_in_PROTAC_literature', 'literature_count_PROTAC', 'pub_id', 'full_id', 'title', 
                                   'Max_halflife', 'Min_halflife',
                                   'Uniprot_keyword', 'Uniprot_PTM', 'Uniprot_CrossLink', 
                                   'Ub_PhosphoSitePlus', 'Ub_mUbiSiDa_2013', 'number_of_ubiquitination_sites', 
                                   'Bcell_mean', 'NKcell_mean', 'Hepatocytes_mean', 'MouseNeuorons_mean',
                                   'count_compound_chembl_ids_PROTAC' #'compound_chembl_ids_PROTAC', #'pchembl_values_PROTAC', 
                                   ]]


        # Score each category, and label highest category
        self.out_df['Clinical_Precedence_PROTAC'] = self.out_df.apply(self._clinical_precedence, axis=1)
        self.out_df['Literature_Precedence_PROTAC'] = 0
        self.out_df.loc[(self.out_df['Bucket_4_PROTAC'] == 1), 'Literature_Precedence_PROTAC'] = 1
        self.out_df['Discovery_Opportunity_PROTAC'] = 0
        self.out_df.loc[((self.out_df['Bucket_5_PROTAC'] == 1) | (self.out_df['Bucket_6_PROTAC'] == 1)) &
                        (self.out_df['Bucket_7_PROTAC'] == 1) & 
                        (self.out_df['Bucket_8_PROTAC'] == 1) & 
                        (self.out_df['PROTAC_location_Bucket'] <= 4),
                        'Discovery_Opportunity_PROTAC'] = 1

        self.out_df['Category_PROTAC'] = 'Unknown'
        self.out_df.loc[((self.out_df['Bucket_5_PROTAC'] == 1) | (self.out_df['Bucket_6_PROTAC'] == 1)) &
                        (self.out_df['Bucket_7_PROTAC'] == 1) & 
                        (self.out_df['Bucket_8_PROTAC'] == 1) & 
                        (self.out_df['PROTAC_location_Bucket'] <= 4),
                        'Category_PROTAC'] = 'Discovery_Opportunity_PROTAC'
        self.out_df.loc[(self.out_df['Bucket_4_PROTAC'] == 1), 'Category_PROTAC'] = 'Literature_Precedence_PROTAC'
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
                                 'Bucket_4_PROTAC': {'literature_count_PROTAC':d.literature_count_PROTAC, 'literature_link':d.full_id, 'title':d.title}, 
                                 'Bucket_5_PROTAC': {'Uniprot_keyword':d.Uniprot_keyword, 'Uniprot_PTM':d.Uniprot_PTM, 'Uniprot_CrossLink':d.Uniprot_CrossLink}, 
                                 'Bucket_6_PROTAC': {'Ub_PhosphoSitePlus':d.Ub_PhosphoSitePlus, 'Ub_mUbiSiDa_2013':d.Ub_mUbiSiDa_2013, 'number_of_ubiquitination_sites':d.number_of_ubiquitination_sites}, 
                                 'Bucket_7_PROTAC': {'Max_halflife':d.Max_halflife, 'Min_halflife':d.Min_halflife}, 
#                                 'Bucket_8_PROTAC': {'compound_chembl_ids_PROTAC':d.compound_chembl_ids_PROTAC}}
                                 'Bucket_8_PROTAC': {'number_chembl_compounds_PROTAC':d.count_compound_chembl_ids_PROTAC}}
            }
