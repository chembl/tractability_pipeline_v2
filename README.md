# Open Targets Tractability Pipeline (version 2)

## Introduction

This pipeline has been developed to produce tractability data for list of input Ensembl Gene IDs. This implementation
is based on the public version of the GSK tractability pipeline, published [here](https://pubs.rsc.org/en/content/articlelanding/2018/md/c7md00633k#!divAbstract)

The pipeline classifies the provided targets into several tractability buckets based on evidence from diverse resources.
Information is provided for different potential modalities, which are subdivided into distinct workflows
(where *X* = bucket number of the specific modality/workflow):
- Small molecule tractability (buckets are denoted with "Bucket_X_sm"),
- Antibody tractability (buckets with "Bucket_X_ab"), 
- PROTAC tractability (buckets with "Bucket_X_PROTAC" and an additional "PROTAC_location_Bucket"), and 
- clinical evidence for other modalities ("Bucket_X_othercl")

The pipeline's main **output** is a TSV file with one target per row, which is additionally converted into a json data structure.



## Installation
Change to the directory containing this file

`pip install .`

Install cxoracle

Set the following environment variables:

`CHEMBL_DB=oracle://address:to@local.chembl` 

`CHEMBL_VERSION=26`

! The current pipeline only works with the most recent ChEMBL version 26 (not older ones) !



## Getting Started

Run the pipeline with the following command:

`run-ot-pipeline-v2 genes.csv`

where `genes.csv` is a file with one Ensembl Gene ID per line with no headers.



## Bucket assignment details

#### Small molecule (SM) buckets
Bucket 1 (ChEMBL) Targets with approved (phase 4) SM drugs  
Bucket 2 (ChEMBL) Targets with SM in clinical phase 2 or 3  
Bucket 3 (ChEMBL) Pre-clinical (phase 1) targets with SM  
Bucket 4 (PDB) Targets with crystal structures with ligands (excluding: organic solvents, sugars & cofactors)  
Bucket 5 (DrugEBIlity) Targets with score ≥0.7  
Bucket 6 (DrugEBIlity) Targets with 0 < score < 0.7  
Bucket 7 (ChEMBL) Targets with ligands (PFI ≤7, SMART hits ≤2, scaffolds ≥2)  
Bucket 8 (Finan et al. 2017) Targets considered druggable using Finan et al's Druggable Genome  
Not implemented: Bucket 9 (SureChEMBL) Targets with ‘chemical’ patents in the last 5 years

#### Antibody (AB) buckets
Bucket 1 (ChEMBL) Targets with approved (phase 4) AB drugs  
Bucket 2 (ChEMBL) Targets with AB in clinical phase 2 or 3  
Bucket 3 (ChEMBL) Pre-clinical (phase 1) targets with AB  
Bucket 4 (Uniprot loc) Targets in “Cell membrane” or “Secreted” - high confidence  
Bucket 5 (GO CC) Targets with GO CC terms indicative for plasma membrane, extracellular region/matrix, or secretion - high confidence  
Bucket 6 (Uniprot loc) Targets in “Cell membrane” or “Secreted” - medium confidence  
Bucket 7 (Uniprot SigP or TMHMM): Targets with predicted Signal Peptide or Trans-membrane regions, and not destined to organelles  
Bucket 8 (GO CC) Targets with GO CC terms indicative for plasma membrane, extracellular region/matrix, or secretion - medium confidence  
Bucket 9 (Human Protein Atlas) Targets with main location “Plasma membrane” - high confidence

#### PROTAC buckets
Bucket 1 (-) Targets with approved (phase 4) PROTAC drugs  
Bucket 2 (-) Targets with PROTAC in clinical phase 2 or 3  
Bucket 3 (-) Pre-clinical (phase 1) targets with PROTAC  
Bucket 4 (Mathieson et al. 2018) Targets with half-life >24 hours  
Bucket 5 (Mathieson et al. 2018) Targets with half-life 10-24 hours  
Bucket 6 (Woong et al. 2011) Targets with known ubiquitination sites  
Not implemented: Bucket 7 (-) Targets with predicted ubiquitination sites  
Bucket 8 (Europe PMC) Targets mentioned in PROTAC literature  
Bucket 9 (ChEMBL) Targets with SM in ChEMBL with activity on target based assay (pChEMBL ≥5) 

PROTAC_location_Bucket  
In addition to PROTAC tractability buckets, there is an additional "PROTAC_location_Bucket", which allows you to assess
whether a target's location is suitable for the PROTAC approach. UniProt and GO locations are assessed. Here, good location = ['cytoplasm', 'cytosol', 'nucleus'], grey location = ['membrane'], and bad location = anywhere else. The PROTAC_location_Bucket score is defined by:  
1 = High confidence good location  
2 = Med confidence good location  
3 = High confidence grey location  
4 = Med condifence grey location  
5 = Unknown location  
6 = Med confidence bad location  
7 = High confidence bad location  

#### Other clinical evidence buckets
Bucket 1 (ChEMBL) Targets with approved (phase4) drugs (other than SM or AB) 
Bucket 2 (ChEMBL) Targets with drug in clinical phase 2 or 3 (other than SM or AB)
Bucket 3 (ChEMBL) Pre-clinical (phase 1) targets with molecule (other than SM or AB)