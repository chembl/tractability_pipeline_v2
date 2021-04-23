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

    def __init__(self, gene_df, store_fetched):

        # list of ensembl IDs for targets to be considered
        self.store_fetched = store_fetched
        self.gene_df = gene_df
        self.gene_list = list(gene_df['ensembl_gene_id'].unique())
        self.human_proteome = self._get_human_proteome()
        self.mapping = self._map_accession_to_gene_id()
        self._merge_genes_mappings_proteome()
        self._add_MyGene_columns()
        
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

    
    def _get_human_proteome(self):
        '''
        Getting the human proteome from UniProt 
        including cross-references to other databases (NCBI GeneID, ChEMBL, BindingDB, DrugBank, PharmGKB, Pharos, PDB) and
        location data (Subcellular location [CC], Transmembrane, Signal peptide, Gene ontology (cellular component));
        Processing entry, gene, and protein names.
        '''
        print("\t- Getting human proteome data from UniProt...")
        #full_url = 'https://www.uniprot.org/uniprot/?query=proteome:UP000005640&format=tab&columns=id,entry%20name,protein%20names,genes(PREFERRED),genes(ALTERNATIVE)'
        full_url = 'https://www.uniprot.org/uniprot/?query=proteome:UP000005640&format=tab&columns=id,entry%20name,protein%20names,genes(PREFERRED),genes(ALTERNATIVE),database(GeneID),database(ChEMBL),database(BindingDB),database(DrugBank),database(PharmGKB),database(Pharos),database(PDB),comment(SUBCELLULAR%20LOCATION),feature(TRANSMEMBRANE),feature(SIGNAL),go(cellular%20component)'

        Uniprot_human_proteome = self.make_request(full_url, data=None)
        Uniprot_human_proteome = [x.split('\t') for x in Uniprot_human_proteome.split('\n')]
        human_proteome = pd.DataFrame(Uniprot_human_proteome[1:], columns=Uniprot_human_proteome[0])
        # rename columns
        human_proteome.rename(columns={'Entry': 'accession',
                                       'Gene names  (primary )': 'gene_name', 
                                       'Gene names  (synonym )': 'gene_name_synonyms'}, inplace=True)
        # only keep row when 'Entry name' is available (discard NAN row)
        human_proteome = human_proteome.loc[human_proteome['Entry name'].notna()]
        # create 'symbol' column
        human_proteome['symbol'] = human_proteome['gene_name']
        #human_proteome[['symbol','Human']] = human_proteome['Entry name'].str.split("_",expand=True)
        # create 'gene_name' column (using first entry in 'Gene names')
        #human_proteome['gene_name'] = human_proteome['Gene names'].str.split(" ",expand=True)[0]
        # create 'protein_name' column (using primary entry, before names in parentheses in 'Protein names', escaping with \\ is reqired)
        human_proteome['protein_name'] = human_proteome['Protein names'].str.split(" \\(",expand=True)[0]
        # save lower case protein name for future case insensitive mapping
        human_proteome['protein_name_lower'] = human_proteome['protein_name'].str.lower()
        # as all protein isoforms have different UniProt IDs, only the first occurence of gene_name is kept 
        # (which should be the primary UniProtID) count: 20487
        human_proteome.drop_duplicates(subset="gene_name", keep='first', inplace=True)
        
        #self.human_proteome = human_proteome
        return human_proteome


    def _map_accession_to_gene_id(self):
        '''
        Mapping UniProt accessions to Ensembl gene IDs (see https://www.uniprot.org/help/api%5Fidmapping)
        '''
        print("\t- Mapping UniProt accessions to Ensembl gene IDs...")
        accession_str = ' '.join(self.human_proteome['accession'].to_list())
        url = 'https://www.uniprot.org/uploadlists/'
        params = {
        'from': 'ACC+ID',
        'to': 'ENSEMBL_ID',
        'format': 'tab',
        'query': accession_str
        }
        
        data = urllib.urlencode(params)
        data = data.encode('utf-8')
        req = urllib2.Request(url, data)
        with urllib2.urlopen(req) as f:
           response = f.read()
        mapping = pd.read_csv(StringIO(response.decode('utf-8')), sep='\t')
        mapping.rename(columns={'From': 'accession', 
                                'To': 'ensembl_gene_id'}, inplace=True)
        return mapping

        
    def _merge_genes_mappings_proteome(self):
        self.id_xref = self.gene_df[['ensembl_gene_id']].merge(self.mapping, how='left')
        self.id_xref = self.id_xref.merge(self.human_proteome, how='left')
        # keep only protein coding genes
        self.id_xref = self.id_xref.loc[~self.id_xref['accession'].isna()]
        
        
    def _add_MyGene_columns(self):
        '''
        Use MyGene to find uniprot accession, entrezgene, pdb, pfam, GO and interpro IDs associated to each Ensembl gene ID.
        :return:
        '''
        print("\t- Getting data from MyGene.Info...")
        mg = mygene.MyGeneInfo()

        # Use MyGene to return list of Uniprot accession numbers for each ensemble gene ID

        mygene_results = mg.getgenes(list(self.gene_list), scopes='ensembl',
                              as_dataframe=True, fields='symbol,uniprot,uniprot.Swiss-Prot,uniprot.TrEMBL,ensembl.protein,entrezgene,pdb,pfam,interpro,interpro.id,go',
                              species='human', returnall=True)

        mygene_results.reset_index(inplace=True)
        mygene_results.rename(columns={'query': 'ensembl_gene_id', 'symbol': 'mygene_symbol'}, inplace=True)
        
        self.id_xref = self.id_xref.merge(mygene_results, how='left', on='ensembl_gene_id')

        if self.store_fetched: 
            self.id_xref.to_csv("{}/id_xref.csv".format(self.store_fetched))

        #self.id_xref.reset_index(inplace=True)
        # self.id_xref.rename({'_id', '_score', 'entrezgene', 'go', 'interpro', 'pdb', 'pfam','uniprot'})
        #self.id_xref['uniprot_accession'] = self.id_xref['uniprot.Swiss-Prot'].apply(self._uniprot_primary_only)
        #self.id_xref.rename(columns={'query': 'ensembl_gene_id'#, 'uniprot_accession': 'accession'
        #                             }, inplace=True)




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
    # drop rows without accession
    out_buckets = out_buckets.dropna(subset=['accession'])
    # drop duplicated accessions
    out_buckets = out_buckets.drop_duplicates(subset="accession", keep='first')

# =============================================================================
#     # =========================================================================
#     # Transforming output DataFrame into json format
#     # =========================================================================
#     if args.workflows_to_run == 'all':
# 
#         with open("./tractability_buckets_{}.json".format(d),"w") as h:
#             for idx,r in out_buckets.iterrows():
#                 gene = {"ensembl_gene_id":r['ensembl_gene_id']}
#                 gene["SM"] = sm.sm2json(r)
#                 gene["AB"] = ab.ab2json(r)
#                 gene["PROTAC"] = protac.protac2json(r)
#                 gene["othercl"] = othercl.othercl2json(r)
#                 json.dump(gene, h)
#     
# 
#         print("\t- as json in './tractability_buckets_{}.json'.\
#                   \n".format(d), sep='\t')
# 
# =============================================================================
   
    # renaming columns
    try:
        out_buckets = out_buckets.rename(columns = {'Bucket_1_sm':'SM_B1_Approved Drug',
                                              'Bucket_2_sm':'SM_B2_Advanced Clinical',
                                              'Bucket_3_sm':'SM_B3_Phase 1 Clinical',
                                              'Bucket_4_sm':'SM_B4_Structure with Ligand',
                                              'Bucket_5_sm':'SM_B5_High-Quality Ligand',
                                              'Bucket_6_sm':'SM_B6_High-Quality Pocket',
                                              'Bucket_7_sm':'SM_B7_Med-Quality Pocket',
                                              'Bucket_8_sm':'SM_B8_Druggable Family',
                                              
                                              'Bucket_1_ab':'AB_B1_Approved Drug',
                                              'Bucket_2_ab':'AB_B2_Advanced Clinical',
                                              'Bucket_3_ab':'AB_B3_Phase 1 Clinical',
                                              'Bucket_4_ab':'AB_B4_UniProt loc high conf',
                                              'Bucket_5_ab':'AB_B5_GO CC high conf',
                                              'Bucket_6_ab':'AB_B6_UniProt loc med conf',
                                              'Bucket_7_ab':'AB_B7_UniProt SigP or TMHMM',
                                              'Bucket_8_ab':'AB_B8_GO CC med conf',
                                              'Bucket_9_ab':'AB_B9_Human Protein Atlas loc',
                                              
                                              'Bucket_1_PROTAC':'PR_B1_Approved Drug',
                                              'Bucket_2_PROTAC':'PR_B2_Advanced Clinical',
                                              'Bucket_3_PROTAC':'PR_B3_Phase 1 Clinical',
                                              'Bucket_4_PROTAC':'PR_B4_Literature',
                                              'Bucket_5_PROTAC':'PR_B5_UniProt Ubiquitination',
                                              'Bucket_6_PROTAC':'PR_B6_Database Ubiquitination',
                                              'Bucket_7_PROTAC':'PR_B7_Half-life Data',
                                              'Bucket_8_PROTAC':'PR_B8_Small Molecule Binder',
                                              'PROTAC_location_Bucket':'PR_locB_Location Score',
                                              
                                              'Bucket_1_othercl':'OC_B1_Approved Drug',
                                              'Bucket_2_othercl':'OC_B2_Advanced Clinical',
                                              'Bucket_3_othercl':'OC_B3_Phase 1 Clinical'
                                             }) #, inplace = True
    except:
        print("Columns could not be renamed.")
    
    # save final output to tsv
    out_buckets.to_csv('tractability_buckets_{}.tsv'.format(d), sep='\t', index=False)
    if store_fetched:
        out_buckets.to_csv('{}/tractability_buckets_{}.tsv'.format(store_fetched,d), sep='\t', index=False)

    # reduced output omitting targets missing all buckets (in all workflows)
    #out_buckets_filled = out_buckets[(out_buckets['Top_bucket_sm'] < 10) | (out_buckets['Top_bucket_ab'] < 10) | (out_buckets['Top_bucket_PROTAC'] < 10) | (out_buckets['Top_bucket_othercl'] < 10)]


    print("\nThe results from all workflows are concatenated and saved\
          \n\t- as table in './tractability_buckets_{}.tsv'".format(d), sep='\t')


if __name__ == '__main__':
    main()

