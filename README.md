# Open Targets Tractability Pipeline (version 2)

## Introduction

This pipeline has been developed to produce tractability data for a list of input Ensembl Gene IDs. This implementation
is based on the public version of the GSK tractability pipeline, published [here](https://pubs.rsc.org/en/content/articlelanding/2018/md/c7md00633k#!divAbstract).

The pipeline assigns the provided targets to tractability buckets based on evidence from diverse resources.
Information is provided for different potential drug modalities, which are subdivided into distinct workflows
(where *X* = bucket number of the specific modality/workflow):
- Small molecule tractability (buckets are denoted with "Bucket_*X*_sm"),
- Antibody tractability (buckets with "Bucket_*X*_ab"), 
- PROTAC tractability (buckets with "Bucket_*X*_PROTAC" and an additional "PROTAC_location_Bucket"), and 
- clinical evidence for other modalities ("Bucket_*X*_othercl")

**Output format**  
The pipeline has two main outputs: 
1) a **TSV** file (table) with one target per row and bucket data as columns, 
2) a **JSON** data structure (generated from the table output) with a hierarchical/nested data structure, relating the evidence to the buckets:  
    {"ensembl_gene_id":{...}, {$modality: {'Bucket_scores':{...}, 'Bucket_evaluation':{...}, 'Category_scores':{...}, 'Category_evaluation':{...}, 'Bucket_evidences':{...}}}}  

Additionaly the data fetched during execution is saved into a 'fetched data' folder (default).


## Installation
Change to the directory containing this file

`pip install .`

Install cxoracle

Set the following environment variables:

`CHEMBL_DB=oracle://address:to@local.chembl` 

`CHEMBL_VERSION=27`

! The current pipeline works with the most recent ChEMBL version 27, and version 26 (but not older ones) !



## Getting Started

Run the pipeline with the following command:

`run-ot-pipeline-v2 genes.csv`

where `genes.csv` is a file with one Ensembl Gene ID per line with no headers.



## Bucket assignment details

#### Small molecule (SM) buckets
Bucket 1 (ChEMBL) Targets with approved (phase 4) SM drug  
Bucket 2 (ChEMBL) Targets with SM in clinical phase 2 or 3  
Bucket 3 (ChEMBL) Targets with SM in clinical phase 1  
Bucket 4 (PDB) Targets with crystal structures with ligands (excluding: buffer compounds, solvents, crystallization agents/additives and sugars)  
Bucket 5 (ChEMBL) Targets with high quality ligands (PFI ≤7, structural_alerts ≤2, pchembl_value ≥5.5 or equivalent activity measure)  
Bucket 6 (DrugEBIlity) Targets with score ≥0.7  
Bucket 7 (DrugEBIlity) Targets with 0 < score < 0.7  
Bucket 8 (Finan et al. 2017) Targets considered druggable using Finan et al's Druggable Genome  
Not implemented: Bucket 9 (SureChEMBL) Targets with ‘chemical’ patents in the last 5 years

#### Antibody (AB) buckets
Bucket 1 (ChEMBL) Targets with approved (phase 4) AB drug  
Bucket 2 (ChEMBL) Targets with AB in clinical phase 2 or 3  
Bucket 3 (ChEMBL) Targets with AB in clinical phase 1  
Bucket 4 (UniProt loc) Targets with UniProt 'Subcellular location' terms indicative for plasma membrane, extracellular region/matrix, or secretion - high confidence  
Bucket 5 (GO CC) Targets with GO 'Cellular component' terms indicative for plasma membrane, extracellular region/matrix, or secretion - high confidence [‘ECO:0000269’ or ‘ECO:0000305’] 
Bucket 6 (UniProt loc) Targets with UniProt 'Subcellular location' terms indicative for plasma membrane, extracellular region/matrix, or secretion - medium confidence  
Bucket 7 (UniProt SigP or TMHMM): Targets with predicted Signal Peptide or Trans-membrane regions, and not destined to organelles  
Bucket 8 (GO CC) Targets with GO 'Cellular component' terms indicative for plasma membrane, extracellular region/matrix, or secretion - medium confidence  
Bucket 9 (Human Protein Atlas) Targets with main location “Plasma membrane” - high confidence

#### PROTAC buckets
Bucket 1 (-) Targets with approved (phase 4) PROTAC drug  
Bucket 2 (-) Targets with PROTAC in clinical phase 2 or 3  
Bucket 3 (-) Targets with PROTAC in clinical phase 1  
Bucket 4 (Europe PMC) Targets mentioned in PROTAC literature (manually curated) + additional information via automated assessments
Bucket 5 (UniProt) Targets with UniProt keyword: "Ubl conjugation [KW-0832]"  
Bucket 6 (PhosphoSitePlus, mUbiSiDa, Woong et al. 2011) Targets with reported ubiquitination sites in PhosphoSitePlus, mUbiSiDa (2013), or [Woong et al. 2011]  
Bucket 7 (Mathieson et al. 2018) Targets with available half-life data  
Bucket 8 (ChEMBL) Targets with SM in ChEMBL with activity on target-based assay (pChEMBL ≥5) 

PROTAC_location_Bucket  
In addition to PROTAC tractability buckets, there is an additional "PROTAC_location_Bucket", which allows you to assess
whether a target's location is suitable for the PROTAC approach. UniProt and GO locations are assessed.  
Here, good location = [cytoplasm or nucleus], grey location = [membrane, including cell membrane and organelle membranes], and bad location = anywhere else.  
The PROTAC_location_Bucket score is defined by:  
1 = High confidence good location  
2 = Med confidence good location  
3 = High confidence grey location  
4 = Med condfidence grey location  
5 = Unknown location  
6 = Med confidence bad location  
7 = High confidence bad location  

#### Other clinical evidence buckets
Bucket 1 (ChEMBL) Targets with approved (phase4) drug (other than SM or AB)  
Bucket 2 (ChEMBL) Targets with drug (other than SM or AB) in clinical phase 2 or 3  
Bucket 3 (ChEMBL) Targets with drug (other than SM or AB) in pre-clinical phase 1  
