#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script runs the whole tractability pipeline (all sub-workflows) and concatenates the final results.
It contains also the functions for the pipeline setup.

Created on Mon Feb  3 10:38:51 2020
@author: Melanie Schneider and Harris Ioannidis
"""

import time
import re
#import zipfile
import zlib
import sys
import argparse
import datetime
import os
import numpy as np
import pandas as pd
import mygene
import pkg_resources
from sqlalchemy import create_engine
import json
import requests
from io import StringIO
import ast
from xml.etree import ElementTree
from urllib.parse import urlparse, parse_qs, urlencode
from requests.adapters import HTTPAdapter, Retry

PY3 = sys.version > '3'
if PY3:
    import urllib.request as urllib2
    import urllib.parse as urllib
else:
    import urllib2
    import urllib

from ot_tractability_pipeline_v2.queries_sm import *
from ot_tractability_pipeline_v2.queries_ab import *
from ot_tractability_pipeline_v2.queries_protac import *
from ot_tractability_pipeline_v2.queries_othercl import *

from ot_tractability_pipeline_v2.buckets_sm import *
from ot_tractability_pipeline_v2.buckets_ab import *
from ot_tractability_pipeline_v2.buckets_protac import *
from ot_tractability_pipeline_v2.buckets_othercl import *


DATA_PATH = pkg_resources.resource_filename('ot_tractability_pipeline_v2', 'data/')

POLLING_INTERVAL = 3
API_URL = "https://rest.uniprot.org"
        
        
retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
session = requests.Session()
session.mount("https://", HTTPAdapter(max_retries=retries))


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
        self._add_IDG_info()

        if self.store_fetched: 
            self.id_xref.to_csv("{}/id_xref.csv".format(self.store_fetched))

        #self.id_xref.reset_index(inplace=True)
        # self.id_xref.rename({'_id', '_score', 'entrezgene', 'go', 'interpro', 'pdb', 'pfam','uniprot'})
        #self.id_xref['uniprot_accession'] = self.id_xref['uniprot.Swiss-Prot'].apply(self._uniprot_primary_only)
        #self.id_xref.rename(columns={'query': 'ensembl_gene_id'#, 'uniprot_accession': 'accession'
        #                             }, inplace=True)
        
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
        # github url now #full_url = 'https://legacy.uniprot.org/uniprot/?query=proteome:UP000005640&format=tab&columns=id,entry%20name,protein%20names,genes(PREFERRED),genes(ALTERNATIVE),database(GeneID),database(ChEMBL),database(BindingDB),database(DrugBank),database(PharmGKB),database(Pharos),database(PDB),comment(SUBCELLULAR%20LOCATION),feature(TRANSMEMBRANE),feature(SIGNAL),go(cellular%20component)'
        #full_url = 'https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Cid%2Cprotein_name%2Cgene_primary%2Cgene_synonym%2Cxref_geneid%2Cxref_chembl%2Cxref_bindingdb%2Cxref_drugbank%2Cxref_pharmgkb%2Cxref_pharos%2Cxref_pdb%2Ccc_subcellular_location%2Cft_transmem%2Cft_signal%2Cgo_c&format=tsv&query=%28proteome:UP000005640%29'
        #full_url = 'https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Cid%2Cprotein_name%2Cgene_primary%2Cgene_synonym%2Cxref_geneid%2Cft_transmem%2Ccc_subcellular_location%2Cgo_c%2Cxref_bindingdb%2Cxref_chembl%2Cxref_drugbank%2Cxref_pdb%2Cxref_pharmgkb%2Cxref_pharos%2Creviewed&format=tsv&query=%28ENSG00000068024%29'
        full_url = 'https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Cid%2Cprotein_name%2Cgene_primary%2Cgene_synonym%2Cxref_geneid%2Cxref_chembl%2Cxref_bindingdb%2Cxref_drugbank%2Cxref_pharmgkb%2Cxref_pharos%2Cxref_pdb%2Ccc_subcellular_location%2Cft_transmem%2Cft_signal%2Cgo_c%2Creviewed&format=tsv&query=%28proteome:UP000005640%29'
        
        Uniprot_human_proteome = self.make_request(full_url, data=None)
        Uniprot_human_proteome = [x.split('\t') for x in Uniprot_human_proteome.split('\n')]
        human_proteome = pd.DataFrame(Uniprot_human_proteome[1:], columns=Uniprot_human_proteome[0])
        human_proteome.query(" `Entry Name` == `Entry Name` & `Reviewed` == 'reviewed' & Entry.str.len() < 7 ", inplace = True)
        human_proteome = human_proteome.drop(columns = ['Reviewed'])
        
        if not human_proteome.empty:
            
            # only keep row when 'Entry name' is available (discard NAN row)
            #human_proteome = human_proteome.loc[human_proteome['Entry Name'].notna()]
        
            # rename columns
            #human_proteome = human_proteome.rename({'Entry': 'accession', 'Gene names  (primary )': 'gene_name', 'Gene names  (synonym )': 'gene_name_synonyms'}, axis = 1)
            human_proteome = human_proteome.rename({'Entry': 'accession', 'Gene Names (primary)': 'gene_name', 'Gene Names (synonym)': 'gene_name_synonyms', 'GeneID': 'Cross-reference (GeneID)', 'ChEMBL': 'Cross-reference (ChEMBL)', 'BindingDB': 'Cross-reference (BindingDB)', 'DrugBank': 'Cross-reference (DrugBank)','PharmGKB': 'Cross-reference (PharmGKB)', 'Pharos': 'Cross-reference (Pharos)', 'PDB':'Cross-reference (PDB)', 'Subcellular location [CC]': 'Subcellular location [CC]', 'Transmembrane':'Transmembrane', 'Gene Ontology (cellular component)': 'Gene ontology (cellular component)'}, axis = 1)
              
            # create 'symbol' column
            human_proteome['symbol'] = human_proteome['gene_name']
            #human_proteome[['symbol','Human']] = human_proteome['Entry name'].str.split("_",expand=True)
            # create 'gene_name' column (using first entry in 'Gene names')
            #human_proteome['gene_name'] = human_proteome['Gene names'].str.split(" ",expand=True)[0]
            # create 'protein_name' column (using primary entry, before names in parentheses in 'Protein names', escaping with \\ is required)
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
        Mapping UniProt accessions to Ensembl gene IDs (see https://www.uniprot.org/help/id_mapping)
        '''
        print("\t- Mapping UniProt accessions to Ensembl gene IDs...")
        human_proteome = self._get_human_proteome()
        ids = human_proteome['accession'].to_list()
        #ids = ' ,'.join(self.human_proteome['accession'].to_list())
        
        n= 100000
        chunks = [ids[i * n:(i + 1) * n] for i in range((len(ids) + n - 1) // n )] 
        
        ##########################################################################################
        #######Functions to be used for mapping Uniprot accession codes to Ensembl IDs############
        ##########################################################################################
        
        def check_response(response):
            try:
                response.raise_for_status()
            except requests.HTTPError:
                print(response.json())
                raise

        def get_next_link(headers):
            re_next_link = re.compile(r'<(.+)>; rel="next"')
            if "Link" in headers:
                match = re_next_link.match(headers["Link"])
                if match:
                    return match.group(1)


        def check_id_mapping_results_ready(job_id):
            while True:
                request = session.get(f"{API_URL}/idmapping/status/{job_id}")
                check_response(request)
                j = request.json()
                if "jobStatus" in j:
                    if j["jobStatus"] == "RUNNING":
                        print(f"Retrying in {POLLING_INTERVAL}s")
                        time.sleep(POLLING_INTERVAL)
                    else:
                        raise Exception(j["jobStatus"])
                else:
                    return bool(j["results"] or j["failedIds"])


        def get_batch(batch_response, file_format, compressed):
            batch_url = get_next_link(batch_response.headers)
            while batch_url:
                batch_response = session.get(batch_url)
                batch_response.raise_for_status()
                yield decode_results(batch_response, file_format, compressed)
                batch_url = get_next_link(batch_response.headers)


        def combine_batches(all_results, batch_results, file_format):
            if file_format == "json":
                for key in ("results", "failedIds"):
                    if key in batch_results and batch_results[key]:
                        all_results[key] += batch_results[key]
            elif file_format == "tsv":
                return all_results + batch_results[1:]
            else:
                return all_results + batch_results
            return all_results


        def get_id_mapping_results_link(job_id):
            url = f"{API_URL}/idmapping/details/{job_id}"
            request = session.get(url)
            check_response(request)
            return request.json()["redirectURL"]


        def decode_results(response, file_format, compressed):
            if compressed:
                decompressed = zlib.decompress(response.content, 16 + zlib.MAX_WBITS)
                if file_format == "json":
                    j = json.loads(decompressed.decode("utf-8"))
                    return j
                elif file_format == "tsv":
                    return [line for line in decompressed.decode("utf-8").split("\n") if line]
                elif file_format == "xlsx":
                    return [decompressed]
                elif file_format == "xml":
                    return [decompressed.decode("utf-8")]
                else:
                    return decompressed.decode("utf-8")
            elif file_format == "json":
                return response.json()
            elif file_format == "tsv":
                return [line for line in response.text.split("\n") if line]
            elif file_format == "xlsx":
                return [response.content]
            elif file_format == "xml":
                return [response.text]
            return response.text


        def get_xml_namespace(element):
            m = re.match(r"\{(.*)\}", element.tag)
            return m.groups()[0] if m else ""


        def merge_xml_results(xml_results):
            merged_root = ElementTree.fromstring(xml_results[0])
            for result in xml_results[1:]:
                root = ElementTree.fromstring(result)
                for child in root.findall("{http://uniprot.org/uniprot}entry"):
                    merged_root.insert(-1, child)
            ElementTree.register_namespace("", get_xml_namespace(merged_root[0]))
            return ElementTree.tostring(merged_root, encoding="utf-8", xml_declaration=True)


        def print_progress_batches(batch_index, size, total):
            n_fetched = min((batch_index + 1) * size, total)
            print(f"Fetched: {n_fetched} / {total}")


        def get_id_mapping_results_search(url):
            parsed = urlparse(url)
            query = parse_qs(parsed.query)
            file_format = query["format"][0] if "format" in query else "json"
            if "size" in query:
                size = int(query["size"][0])
            else:
                size = 500
                query["size"] = size
            compressed = (
                query["compressed"][0].lower() == "true" if "compressed" in query else False
            )
            parsed = parsed._replace(query=urlencode(query, doseq=True))
            url = parsed.geturl()
            request = session.get(url)
            check_response(request)
            results = decode_results(request, file_format, compressed)
            total = int(request.headers["x-total-results"])
            print_progress_batches(0, size, total)
            for i, batch in enumerate(get_batch(request, file_format, compressed), 1):
                results = combine_batches(results, batch, file_format)
                print_progress_batches(i, size, total)
            if file_format == "xml":
                return merge_xml_results(results)
            return results


        def get_id_mapping_results_stream(url):
            if "/stream/" not in url:
                url = url.replace("/results/", "/results/stream/")
            request = session.get(url)
            check_response(request)
            parsed = urlparse(url)
            query = parse_qs(parsed.query)
            file_format = query["format"][0] if "format" in query else "json"
            compressed = (
                query["compressed"][0].lower() == "true" if "compressed" in query else False
            )
            return decode_results(request, file_format, compressed)
        
        def submit_id_mapping(from_db, to_db, ids):
            request = requests.post(
                f"{API_URL}/idmapping/run",
                data={"from": from_db, "to": to_db, "ids": ",".join(ids)},
            )
            check_response(request)
            return request.json()["jobId"]
        
        ##########################################################################################
        
        
        mappings_ = []
        
        for chk in chunks:
            
            job_id = submit_id_mapping(from_db = "UniProtKB_AC-ID", to_db = "Ensembl", ids = chk)

            if check_id_mapping_results_ready(job_id):
                link = get_id_mapping_results_link(job_id)
                results = get_id_mapping_results_search(link)
                
            for item in results.items():
        
                mapping_ = pd.DataFrame.from_dict(item[1])
                mapping_.rename(columns = {'from': 'accession', 'to': 'ensembl_gene_id'}, inplace = True)
                mappings_.append(mapping_)
            
        mapping_ = pd.concat(mappings_)
        mapping_
        
        if not mapping_.empty:
            
           # mapping_.rename(columns={'From': 'accession', 'To': 'ensembl_gene_id'}, inplace=True)
            mapping_['ensembl_gene_id'] = mapping_['ensembl_gene_id'].str.split("." , expand = True)[0]
            mapping = mapping_

        return mapping

        
    def _merge_genes_mappings_proteome(self):
        self.id_xref = self.gene_df[['ensembl_gene_id']].merge(self.mapping, how='left')
        self.id_xref = self.id_xref.merge(self.human_proteome, how='left')
        # keep only protein coding genes
        self.id_xref = self.id_xref.loc[~self.id_xref['accession'].isna()]
        
        
    @staticmethod
    def extract_terms(dfcolumn):
        '''
        Extract GO terms from column returned by MyGene query.
        '''
        go_terms_list = []
        try:
            terms = dfcolumn
            #terms = ast.literal_eval(terms)
        except:
            return []
        if isinstance(terms, float):
            return []
        if isinstance(terms, dict):
            terms = [terms]
        for terms_dict in terms:
            try:
                go_terms = terms_dict['term']
                go_terms_list.append(go_terms)
            except TypeError:
                continue
        return go_terms_list

        
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

        # check if columns exist and if not create them (filled with np.NaN or "")
        cols1 = ['go.BP','go.CC','go.MF']
        cols2 = ['go.BP.evidence','go.CC.evidence','go.MF.evidence','go.BP.gocategory','go.CC.gocategory','go.MF.category','go.BP.id','go.CC.id','go.MF.id','go.BP.qualifier','go.CC.qualifier','go.MF.qualifier','go.BP.term','go.CC.term','go.MF.term']
        self.id_xref = self.id_xref.reindex(self.id_xref.columns.union(cols1, sort=False), axis=1, fill_value=np.NaN)
        self.id_xref = self.id_xref.reindex(self.id_xref.columns.union(cols2, sort=False), axis=1, fill_value="")
        # fix missing data in aggregated go columns (When there is only one entry for BP/CC/MF, the aggregated go.BP/CC/MF column is NAN, but the data is available in the seperate respective columns and thus fetched from there.)
        self.id_xref['go.BP'].fillna("[{'evidence': '"+self.id_xref['go.BP.evidence']+"', 'gocategory': '"+self.id_xref['go.BP.gocategory']+"', 'id': '"+self.id_xref['go.BP.id']+"', 'qualifier': '"+self.id_xref['go.BP.qualifier']+"', 'term': '"+self.id_xref['go.BP.term']+"'}]", inplace=True)
        self.id_xref['go.CC'].fillna("[{'evidence': '"+self.id_xref['go.CC.evidence']+"', 'gocategory': '"+self.id_xref['go.CC.gocategory']+"', 'id': '"+self.id_xref['go.CC.id']+"', 'qualifier': '"+self.id_xref['go.CC.qualifier']+"', 'term': '"+self.id_xref['go.CC.term']+"'}]", inplace=True)
        self.id_xref['go.MF'].fillna("[{'category': '"+self.id_xref['go.MF.category']+"', 'evidence': '"+self.id_xref['go.MF.evidence']+"', 'id': '"+self.id_xref['go.MF.id']+"', 'qualifier': '"+self.id_xref['go.MF.qualifier']+"', 'term': '"+self.id_xref['go.MF.term']+"'}]", inplace=True)
        
        self.id_xref['GO_BioProcess'] = self.id_xref['go.BP'].apply(self.extract_terms)
        self.id_xref['GO_MolFunction'] = self.id_xref['go.MF'].apply(self.extract_terms)
        
        
    def _add_IDG_info(self):
        '''
        Add information from Illuminating the Druggable Genome (IDG) project - Pharos/TCRD: protein family classification and Target Development Levels (TDLs)
        '''
        print("\t- Getting data from IDG/TCRD version 6.11.0...")
        conn_string = 'mysql://tcrd@tcrd.kmc.io:3306/tcrd6110'
        # For latest database version check http://juniper.health.unm.edu/tcrd/
        # db_connection = sql.connect(host='tcrd.kmc.io', db='tcrd6110', user='tcrd')
        engine = create_engine(conn_string)
        # Read in everything from the protein table
        query1 = "SELECT id, name, description, uniprot, sym, family ,dtoclass \
                 FROM protein"
        protein = pd.read_sql(query1, con=engine)
        # Read in info from the target table
        query2 = "SELECT id, name, ttype, tdl, fam \
                 FROM target" # \ WHERE tdl='Tclin' OR tdl='Tchem'"
        target = pd.read_sql(query2, con=engine)
        # Closing the connection
        engine.dispose()
        # Joining the Target and Protein Data Frames on internal id
        protein = protein.set_index("id")
        target = target.set_index("id")
        # Concatenate protein and target dataframes (excluding column 'name' from target):
        result = pd.concat([protein, target.drop('name', axis=1)], axis=1, join='outer')
        result.drop_duplicates(subset="uniprot", keep='first', inplace=True)
        # added columns are: name, description, uniprot, sym	, family	, dtoclass, 	ttype, tdl, 	fam
        result.rename(columns={'family': 'IDG_family', 
                               'dtoclass': 'IDG_dtoclass', 
                               'ttype': 'IDG_ttype', 
                               'tdl': 'IDG_tdl', 
                               'fam': 'IDG_fam'}, inplace=True)
        # Merge IDG results with self.id_xref
        self.id_xref = self.id_xref.merge(result, how='left', left_on='accession', right_on='uniprot')




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
    gene_df = pd.read_csv(args.genes, encoding='utf-8', header=None, names=['ensembl_gene_id'], usecols=[0])
    gene_df = gene_df.drop_duplicates(subset=['ensembl_gene_id'])
    
    # URL to local ChEMBL database
    global database_url
    database_url = args.db
    # If database url not supplied, get from environemnt variable
    if database_url is None:
        database_url = os.getenv('CHEMBL_DB')
    print("Using ChEMBL database at {}".format(database_url))

    # check workflows_to_run 
    if args.workflows_to_run not in ['all','sm_ab_othercl','sm_protac_othercl','sm_othercl']:
        parser.print_help(sys.stderr)
        sys.exit(1)
    else: print("Workflows to be run are: {}".format(args.workflows_to_run))
    
    # create external data folder
    if args.store_fetched_data == True:
        ts = datetime.datetime.now().strftime("%Y-%m-%d %Hh%Mm%Ss")
        fn = "./fetched data {}".format(ts)
        os.mkdir(fn)
        store_fetched = fn
        print("Additional data will be saved to {}".format(store_fetched))
    else: store_fetched = False 
    
# =============================================================================
#     run(ensembl_id_list, database_url, out_file_name=args.out_file, store_fetched=store_fetched)
# 
# def run(ensembl_id_list, database_url, out_file_name, store_fetched):
# =============================================================================    

    # Assign tractability buckets
    
    d=datetime.date.today()

    # From Ensembl gene ID, get Uniprot, GO.CC, Symbol and PDB codes
    print("\nSetting up the pipeline...\n")
    setup = Pipeline_setup(gene_df, store_fetched)

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

