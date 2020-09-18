#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script runs the whole tractability pipeline (all sub-workflows) and concatenates the final results.
It contains also the functions for the pipeline setup.

Created on Mon Feb  3 10:38:51 2020
@author: Melanie Schneider
"""

#import time
#import zipfile
import sys
import argparse
import datetime
import os

import mygene
import numpy as np
import pandas as pd
import pkg_resources
from sqlalchemy import create_engine
import json
import requests

PY3 = sys.version > '3'
if PY3:
    import urllib.request as urllib2
else:
    import urllib2

from ot_tractability_pipeline_v2.queries_sm import *
from ot_tractability_pipeline_v2.queries_ab import *
from ot_tractability_pipeline_v2.queries_protac import *
from ot_tractability_pipeline_v2.queries_othercl import *

from ot_tractability_pipeline_v2.buckets_sm import *
from ot_tractability_pipeline_v2.buckets_ab import *
from ot_tractability_pipeline_v2.buckets_protac import *
from ot_tractability_pipeline_v2.buckets_othercl import *



DATA_PATH = pkg_resources.resource_filename('ot_tractability_pipeline_v2', 'data/')



class Pipeline_setup(object):
    '''
    Information retrieved from MyGene is used in both the small molecule and antibody pipelines. This class handles
    the aggregation of data ahead of assigning targets to buckets
    '''

    def __init__(self, ensembl_gene_id_list, store_fetched):

        # list of ensembl IDs for targets to be considered
        self.store_fetched = store_fetched
        self.gene_list = ensembl_gene_id_list
        self._add_MyGene_columns()
        self._add_biomart_columns()
#        self._add_ebi_columns()
        self._accession_column()
        
        # # If database url not supplied, get from environemnt variable
        # if database_url is None:
        #     database_url = os.getenv('CHEMBL_DB')

        # Create ChEMBL DB connection
        try:
            self.engine = create_engine(database_url)
        except AttributeError:
            print('''
            Please supply a valid database URL to your local ChEMBL installation using one of the following methods
            
            1) Supply the URL to this class within your script
            2) Create an environment variable named 'CHEMBL_DB' with containing your database URL
            3) Use the command line flag '--db' if using the supplied run script  '''
                  )
            raise


    def _uniprot_primary_only(self, s):
        '''
        If multiple uniprot IDs, only take primary (assume first in list is primary)
        :return:
        '''

        if isinstance(s, list):
            return s[0]
        else:
            return s

    def _add_MyGene_columns(self):
        '''
        Use MyGene to find uniprot accession, entrezgene  pdb, pfam, GO and interpro IDs associated to each Ensembl gene ID.
        :return:
        '''

        mg = mygene.MyGeneInfo()

        # Use MyGene to return list of Uniprot accession numbers for each ensemble gene ID

        results = mg.getgenes(list(self.gene_list), scopes='ensembl',
                              as_dataframe=True, fields='symbol,uniprot,uniprot.Swiss-Prot,uniprot.TrEMBL,ensembl.protein,entrezgene,pdb,pfam,interpro,interpro.id,go',
                              species='human', returnall=True)

        self.id_xref = results

        self.id_xref['uniprot_accession'] = self.id_xref['uniprot.Swiss-Prot'].apply(self._uniprot_primary_only)
        self.id_xref.reset_index(inplace=True)
        # self.id_xref.rename({'_id', '_score', 'entrezgene', 'go', 'interpro', 'pdb', 'pfam','uniprot'})
        self.id_xref.rename(columns={'query': 'ensembl_gene_id'#, 'uniprot_accession': 'accession'
                                     }, inplace=True)


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

    
    def _add_biomart_columns(self):
        '''
        Use BioMart to find uniprot accession IDs associated to each Ensembl gene ID. 
        For more info see: https://www.ensembl.org/info/data/biomart/biomart_restful.html
        The query is in xml format:
            <?xml version="1.0" encoding="UTF-8"?>
            <!DOCTYPE Query>
            <Query  virtualSchemaName = "default" formatter = "TSV" header = "1" uniqueRows = "1" count = "" datasetConfigVersion = "0.6" >
                <Dataset name = "hsapiens_gene_ensembl" interface = "default" >
                    <Attribute name = "ensembl_gene_id" />
                    <Attribute name = "hgnc_symbol" />
                    <Attribute name = "uniprotswissprot" />
                </Dataset>
            </Query>
        
        :return:
        '''

        full_url = 'http://www.ensembl.org/biomart/martservice?query=%3C?xml%20version=%221.0%22%20encoding=%22UTF-8%22?%3E%3C!DOCTYPE%20Query%3E%3CQuery%20%20virtualSchemaName%20=%20%22default%22%20formatter%20=%20%22TSV%22%20header%20=%20%221%22%20uniqueRows%20=%20%220%22%20count%20=%20%22%22%20datasetConfigVersion%20=%20%220.6%22%20%3E%3CDataset%20name%20=%20%22hsapiens_gene_ensembl%22%20interface%20=%20%22default%22%20%3E%3CAttribute%20name%20=%20%22ensembl_gene_id%22%20/%3E%3CAttribute%20name%20=%20%22hgnc_symbol%22%20/%3E%3CAttribute%20name%20=%20%22uniprotswissprot%22%20/%3E%3C/Dataset%3E%3C/Query%3E'
        
        #from ot_tractability_pipeline_v2.buckets_ab import Antibody_buckets        #biomart = Antibody_buckets.make_request(full_url, data=None)
        biomart = self.make_request(full_url, data=None)
        biomart = [x.split('\t') for x in biomart.split('\n')]
        biomart = pd.DataFrame(biomart[1:], columns=biomart[0])
        biomart.rename(columns={'Gene stable ID': 'ensembl_gene_id', 
                                'HGNC symbol': 'hgnc_symbol', 
                                'UniProtKB/Swiss-Prot ID': 'uniprotswissprot'}, inplace=True)
        # replace empty string with NaN
        biomart = biomart.replace(r'^\s*$', np.NaN, regex=True)
        # only keep row when 'uniprotswissprot' is available (discard NAN row)
        biomart = biomart.loc[biomart['uniprotswissprot'].notna()]
        
        self.id_xref = biomart.merge(self.id_xref, how='right', on='ensembl_gene_id')
        self.id_xref.drop_duplicates(['ensembl_gene_id'], inplace=True, ignore_index=True)        

    
    def _get_IDs_from_json(self, item):
        return {'ensembl_gene_id': item['uniProtKBCrossReferences'][0]['properties'][1]['value'],
                'ebi_accession': item['primaryAccession'],
                'isoforms': "Yes" if 'isoformId' in item['uniProtKBCrossReferences'][0] else "No",
                #'ensembl_prot_id': item['uniProtKBCrossReferences'][0]['properties'][0]['value'],
                #'isoformId': item['uniProtKBCrossReferences'][0].get('isoformId', None),
                }

# =============================================================================
#     def _add_ebi_columns(self):
#         '''
#         Use EBI UniProt dev API to find most recent/up-to-date uniprot accession IDs associated to each Ensembl gene ID.
#         :return:
#         '''
#         
#         requestURL = "https://wwwdev.ebi.ac.uk/uniprot/api/uniprotkb/download?fields=xref_ensembl&query=(database:ensembl)%20AND%20(organism_id:9606)&size=100000"
#         response = requests.get(requestURL, headers={"Accept" : "application/json"})
#         
#         data = response.json()
#         data = data['results']
#         ebi_df = pd.DataFrame([self._get_IDs_from_json(item) for item in data])
#         
#         # keep only first (mostly primary) UniProt ID
#         ebi_df_unique = ebi_df.drop_duplicates(['ensembl_gene_id'])
#         
#         self.id_xref = ebi_df_unique.merge(self.id_xref, how='right', on='ensembl_gene_id')
#         self.id_xref.drop_duplicates(['ensembl_gene_id'], inplace=True, ignore_index=True)
# 
# =============================================================================

    def _accession_column(self):
        '''
        Create accession colum with values from MyGene uniprot_accession and replace missing entries 
        with BioMart uniprotswissprot, and if still empty with ebi_accession 
        (ebi_accession has more annotations than MyGene, but it is not always the primary one)
        :return:
        '''
        self.id_xref.insert(loc=1, column='accession', value=self.id_xref['uniprot_accession'].fillna(value=self.id_xref['uniprotswissprot']))#.fillna(value=self.id_xref['ebi_accession']))
        
        print(self.id_xref.columns)
        
        if self.store_fetched: 
            self.id_xref.to_csv("{}/id_xref.csv".format(self.store_fetched))



#--------------------MAIN-----------------------
            
def main(args=None):
    
    args = None
    if args is None:
        #--------------------------------
        ## for debugging:
        #sys.argv = 'script.py /Users/melanie/tractability_project/tractability_pipeline_v2/input_genelists/20.04_target_list.csv'.split()
        #sys.argv = 'script.py /Users/melanie/tractability_project/tractability_pipeline_v2/input_genelists/test01.csv'.split()
        #--------------------------------
        args = sys.argv[1:]
    
    
    parser = argparse.ArgumentParser(description='Assess tractability of gene list')
    parser.add_argument('genes',
                        help='A file of Ensembl Gene IDs, one per line and no header')
    parser.add_argument('--workflows_to_run', default='all',
                        help='Workflows to be run - options are: "all"(default), "sm_ab_othercl", "sm_protac_othercl" and "sm_othercl"')
    parser.add_argument('--db', default=None,
                        help='Address to your local ChEMBL installation. Preferred method is to set the CHEMBL_DB environment variable')    
    parser.add_argument('--store_fetched_data', default=True,
                        help='Store data fetched from external resources, good habit for future troubleshooting - bool, default True')
#    parser.add_argument('--out_file', default='tractability_buckets.tsv',
#                        help='Name of output csv file, default "tractability_buckets.tsv"')
    
    args = parser.parse_args()
    
    
    
    
    # Get a unique list of Gene IDs from OT
    gene_list = pd.read_csv(args.genes, encoding='utf-8', header=None, names=['ensembl_gene_id'], usecols=[0])
    ensembl_id_list = list(gene_list['ensembl_gene_id'].unique())
    
    # URL to local ChEMBL database
    global database_url
    database_url = args.db
    # If database url not supplied, get from environemnt variable
    if database_url is None:
        database_url = os.getenv('CHEMBL_DB')
    
    # create external data folder
    if args.store_fetched_data == True:
        ts = datetime.datetime.now().strftime("%Y-%m-%d %Hh%Mm%Ss")
        fn = "./fetched data {}".format(ts)
        os.mkdir(fn)
        store_fetched = fn
    else: store_fetched = False 

    # check workflows_to_run 
    if args.workflows_to_run not in ['all','sm_ab_othercl','sm_protac_othercl','sm_othercl']:
        parser.print_help(sys.stderr)
        sys.exit(1)
    
# =============================================================================
#     run(ensembl_id_list, database_url, out_file_name=args.out_file, store_fetched=store_fetched)
# 
# def run(ensembl_id_list, database_url, out_file_name, store_fetched):
# =============================================================================
    

    # Assign tractability buckets
    
    d=datetime.date.today()

    # From Ensembl gene ID, get Uniprot, GO.CC, Symbol and PDB codes
    print("\nSetting up the pipeline...\n")
    setup = Pipeline_setup(ensembl_id_list, store_fetched)

    # Get small molecule tractability info
    print("\nRunning the Small Molecule workflow...\n")
    sm = Small_molecule_buckets(setup, database_url=database_url)
    out_buckets_sm = sm.assign_buckets()
    print(out_buckets_sm.groupby('Top_bucket_sm')['accession'].count())
    if store_fetched: out_buckets_sm.to_csv('{}/tractability_buckets_{}_sm.tsv'.format(store_fetched,d), sep='\t')
    print("\nThe Small Molecule workflow run successfully. \nThe results are saved in '{}/tractability_buckets_{}_sm.tsv'.\n".format(store_fetched,d), sep='\t')

    # Get antibody tractability info
    if args.workflows_to_run in ['all','sm_ab_othercl']:
        print("\nRunning the Antibody workflow...\n")
        ab = Antibody_buckets(setup, database_url=database_url)
        out_buckets_ab = ab.assign_buckets()
        print(out_buckets_ab.groupby('Top_bucket_ab')['accession'].count())
        if store_fetched: out_buckets_ab.to_csv('{}/tractability_buckets_{}_ab.tsv'.format(store_fetched,d), sep='\t')
        print("\nThe Antibody workflow run successfully. \nThe results are saved in '{}/tractability_buckets_{}_ab.tsv'.\n".format(store_fetched,d), sep='\t')

    # Get protac tractability info
    if args.workflows_to_run in ['all','sm_protac_othercl']:
        print("\nRunning the PROTAC workflow...\n")
        protac = Protac_buckets(setup, database_url=database_url)
        out_buckets_protac = protac.assign_buckets()
        print(out_buckets_protac.groupby('Top_bucket_PROTAC')['accession'].count())
        if store_fetched: out_buckets_protac.to_csv('{}/tractability_buckets_{}_protac.tsv'.format(store_fetched,d), sep='\t')
        print("\nThe PROTAC workflow run successfully. \nThe results are saved in '{}/tractability_buckets_{}_protac.tsv'.\n".format(store_fetched,d), sep='\t')

    # Get other clinical info
    if args.workflows_to_run in ['all','sm_ab_othercl','sm_protac_othercl','sm_othercl']:
        print("\nRunning the OtherClinical workflow...\n")
        othercl = Othercl_buckets(setup, database_url=database_url)
        out_buckets_othercl = othercl.assign_buckets()
        print(out_buckets_othercl.groupby('Top_bucket_othercl')['accession'].count())
        if store_fetched: out_buckets_othercl.to_csv('{}/tractability_buckets_{}_othercl.tsv'.format(store_fetched,d), sep='\t')
        print("\nThe OtherClinical workflow run successfully. \nThe results are saved in '{}/tractability_buckets_{}_othercl.tsv'.\n".format(store_fetched,d), sep='\t')



    
    # =============================================================================
    # Combine output from different workflows to one final summary output
    # =============================================================================
    # try:
    # if out_buckets_sm and out_buckets_ab and out_buckets_protac and out_buckets_othercl:
    
    #data_frames = [out_buckets_sm, out_buckets_ab, out_buckets_protac, out_buckets_othercl]
    if args.workflows_to_run == 'all':
        other_frames = [out_buckets_ab.drop(columns=['accession', 'symbol']), out_buckets_protac.drop(columns=['accession', 'symbol']), out_buckets_othercl.drop(columns=['accession', 'symbol'])]
#        data_frames = [out_buckets_sm, out_buckets_ab, out_buckets_protac, out_buckets_othercl]
    elif args.workflows_to_run == 'sm_ab_othercl':
        other_frames = [out_buckets_ab.drop(columns=['accession', 'symbol']), out_buckets_othercl.drop(columns=['accession', 'symbol'])]
#        data_frames = [out_buckets_sm, out_buckets_ab, out_buckets_othercl]
    elif args.workflows_to_run == 'sm_protac_othercl':
        other_frames = [out_buckets_protac.drop(columns=['accession', 'symbol']), out_buckets_othercl.drop(columns=['accession', 'symbol'])]
#        data_frames = [out_buckets_sm, out_buckets_protac, out_buckets_othercl]
    elif args.workflows_to_run == 'sm_othercl':
        other_frames = [out_buckets_othercl.drop(columns=['accession', 'symbol'])]
#        data_frames = [out_buckets_sm, out_buckets_othercl]
    
    # combining workflow output data frames
#    out_buckets = pd.concat(data_frames, axis=1, join='outer', sort=False)

    out_buckets = out_buckets_sm.join(other_frames, how='left', sort=False)

    out_buckets.reset_index(inplace=True)
    out_buckets.rename(columns={'index': 'ensembl_gene_id'}, inplace=True)
    
    # Remove duplicated columns (such as 'symbol', 'accession')
    out_buckets = out_buckets.loc[:,~out_buckets.columns.duplicated()]

    # sorting output dataframe
    if args.workflows_to_run == 'all':
        out_buckets.sort_values(['Clinical_Precedence_sm', 'Discovery_Precedence_sm', 'Predicted_Tractable_sm',
                                  'Clinical_Precedence_ab', 'Predicted_Tractable_ab_High_confidence', 'Predicted_Tractable_ab_Medium_to_low_confidence',
                                  'Top_bucket_PROTAC', 'PROTAC_location_Bucket',
                                  'Clinical_Precedence_othercl'],
                                ascending=[False, False, False,
                                           False, False, False,
                                           True, True, 
                                           False], inplace=True)

    if args.workflows_to_run == 'sm_ab_othercl':
        out_buckets.sort_values(['Clinical_Precedence_sm', 'Discovery_Precedence_sm', 'Predicted_Tractable_sm',
                                  'Clinical_Precedence_ab', 'Predicted_Tractable_ab_High_confidence', 'Predicted_Tractable_ab_Medium_to_low_confidence',
                                  'Clinical_Precedence_othercl'],
                                ascending=[False, False, False,
                                           False, False, False,
                                           False], inplace=True)

    if args.workflows_to_run == 'sm_protac_othercl':
        out_buckets.sort_values(['Clinical_Precedence_sm', 'Discovery_Precedence_sm', 'Predicted_Tractable_sm',
                                  'Top_bucket_PROTAC', 'PROTAC_location_Bucket',
                                  'Clinical_Precedence_othercl'],
                                ascending=[False, False, False,
                                           True, True, 
                                           False], inplace=True)

    if args.workflows_to_run == 'sm_othercl':
        out_buckets.sort_values(['Clinical_Precedence_sm', 'Discovery_Precedence_sm', 'Predicted_Tractable_sm',
                                  'Clinical_Precedence_othercl'],
                                ascending=[False, False, False,
                                           False], inplace=True)

    # drop rows without gene symbol
    out_buckets = out_buckets.dropna(subset=['symbol'])
    # drop duplicated accessions
    out_buckets = out_buckets.drop_duplicates(subset="accession", keep='first')
    
    # save final output to tsv
    out_buckets.to_csv('tractability_buckets_{}.tsv'.format(d), sep='\t', index=False)
    if store_fetched:
        out_buckets.to_csv('{}/tractability_buckets_{}.tsv'.format(store_fetched,d), sep='\t', index=False)

    # drop entries without UniProt 'accession'
    out_buckets = out_buckets.dropna(subset=['accession']) #, inplace=True

    # reduced output omitting targets missing all buckets (in all workflows)
    #out_buckets_filled = out_buckets[(out_buckets['Top_bucket_sm'] < 10) | (out_buckets['Top_bucket_ab'] < 10) | (out_buckets['Top_bucket_PROTAC'] < 10) | (out_buckets['Top_bucket_othercl'] < 10)]


    print("\nThe results from all workflows are concatenated and saved\
          \n\t- as table in './tractability_buckets_{}.tsv'".format(d), sep='\t')

    # =========================================================================
    # Transforming output DataFrame into json format
    # =========================================================================
    if args.workflows_to_run == 'all':

        with open("./tractability_buckets_{}.json".format(d),"w") as h:
            for idx,r in out_buckets.iterrows():
                gene = {"ensembl_gene_id":r['ensembl_gene_id']}
                gene["SM"] = sm.sm2json(r)
                gene["AB"] = ab.ab2json(r)
                gene["PROTAC"] = protac.protac2json(r)
                gene["othercl"] = othercl.othercl2json(r)
                json.dump(gene, h)
    

        print("\t- as json in './tractability_buckets_{}.json'.\
                  \n".format(d), sep='\t')

                 


if __name__ == '__main__':
    main()

