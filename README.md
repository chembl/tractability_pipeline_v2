# Open Targets Tractability Pipeline (version 2)

## Introduction

This pipeline has been developed to produce tractability data for list of input Ensembl Gene IDs. This implementation
is based on the public version of the GSK tractability pipeline, published [here](https://pubs.rsc.org/en/content/articlelanding/2018/md/c7md00633k#!divAbstract)

The pipeline classifies the provided targets into several tractability buckets based on evidence from diverse resources.
Information is provided for different potential modalities
(where *X* = bucket number of the specific modality):
- Small molecule tractability (buckets are denoted with "Bucket_X_sm"),
- Antibody tractability (buckets with "Bucket_X_ab"), 
- PROTAC tractability (buckets with "Bucket_X_PROTAC" and an additional "PROTAC_location_Bucket"), and 
- clinical evidence for other modalities ("Bucket_X_other")

The pipeline's main **output** is a TSV file with one target per row. 



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

#### Small molecule buckets
Bucket 1 (ChEMBL) Targets with approved SM drugs (phase4)  
Bucket 2 (ChEMBL) Targets with SM in ≥phase2  
Bucket 3 (ChEMBL) Pre-clinical targets with SM  
Bucket 4 (PDB) Targets with crystal structures with ligands (excluding: organic solvents, sugars & cofactors)  
Bucket 5 (DrugEBIlity) Targets with score ≥0.7  
Bucket 6 (DrugEBIlity) Targets with 0 < score < 0.7  
Bucket 7 (ChEMBL) Targets with ligands (PFI ≤7, SMART hits ≤2, scaffolds ≥2)  
Bucket 8 (Finan et al. 2017) Targets considered druggable using Finan et al's druggable genome  
Bucket 9 (SureChEMBL) Targets with ‘chemical’ patents in the last 5 years  

#### Antibody buckets
Bucket 1 (ChEMBL) Targets with approved AB drugs (phase4)  
Bucket 2 (ChEMBL) Targets with AB in ≥phase2  
Bucket 3 (ChEMBL) Pre-clinical targets with AB  
Bucket 4 (HPA) Targets in “Plasma membrane”, high confidence  
Bucket 5 (Uniprot loc) Targets in “Cell membrane” or “Secreted”, high confidence  
Bucket 6 (Uniprot loc) Targets in “Cell membrane” or “Secreted” or “Membrane”, low or unknown confidence  
Bucket 7 (SigP + TMHMM) Targets with predicted signal peptide or trans-membrane regions, and not destined to organelles  
Bucket 8 (GO CC) Targets with the parent term GO: 0005576 (extracellular region) or GO: 0031012 (extracellular matrix) or GO: 0005886 (plasma membrane) or their child terms.  

#### PROTAC buckets
Bucket 1 (-) Targets with approved PROTAC drugs (phase4)  
Bucket 2 (-) Targets with PROTAC in ≥phase2  
Bucket 3 (-) Pre-clinical targets with PROTAC  
Bucket 4 (Mathieson et al. 2018) Targets with half-life >100 hours  
Bucket 5 (Mathieson et al. 2018) Targets with half-life 24-100 hours  
Bucket 6 (Woong et al. 2011) Targets with known ubiquitination sites  
Bucket 7 (-) Targets with predicted ubiquitination sites  
Bucket 8 (Europe PMC) Targets mentioned in PROTAC literature  
Bucket 9 (ChEMBL) Target with SM in ChEMBL with activity on target based assay  

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
