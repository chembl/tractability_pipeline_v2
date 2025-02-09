# Open Targets Tractability Pipeline (version 2)

## Introduction

This tractability pipeline has been developed to produce tractability data for a list of input Ensembl Gene IDs. This implementation
is based on the public version of the GSK tractability pipeline, published [here](https://pubs.rsc.org/en/content/articlelanding/2018/md/c7md00633k#!divAbstract), 
which assesses proteins of interest against a series of criteria leading to them being assigned to a tractability bucket based on either literature or derived knowledge. 
This formalism formed the basis of the tractability schema recently implemented in the Open Targets platform.

The Open Targets tractability pipeline assigns the provided targets to tractability buckets (usually 1 or 0) based on evidence from diverse resources.
Information is provided for different potential drug modalities, which are subdivided into distinct workflows
(where *X* = bucket number of the specific modality workflow):
- Small molecule tractability (buckets are denoted with "Bucket_*X*_sm"),
- Antibody tractability (buckets with "Bucket_*X*_ab"), 
- PROTAC tractability (buckets with "Bucket_*X*_PROTAC" and an additional "PROTAC_location_Bucket"), and 
- clinical evidence for other modalities ("Bucket_*X*_othercl")

**Output format**  
The pipeline has one main output:
- a **TSV** file (table) with one target per row and bucket data as columns.

The column names of tractability buckets are replaced with more descriptive names prior to saving the output file:

  'Bucket_1_sm':'SM_B1_Approved Drug',  
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


It is possible to activate a second output format by uncommenting the respective lines in run_pipeline.py):  
- a **JSON** data structure (generated from the table output) with a hierarchical/nested data structure, relating the evidence to the buckets:  
    {"ensembl_gene_id":{...}, {$modality: {'Bucket_scores':{...}, 'Bucket_evaluation':{...}, 'Category_scores':{...}, 'Category_evaluation':{...}, 'Bucket_evidences':{...}}}}  

Additionaly the data fetched during execution is saved into a 'fetched data' folder (default).


## Installation

To install the pipeline, change to the directory containing this README file (contains also the setup.py file), then run:

`pip install .`

From the same directory, install the PROTAC named entity recognition (NER) model by running:

`pip install ./ot_tractability_pipeline_v2/en_NER_PROTAC-0.2.5.tar.gz`

Install cx-oracle (if an Oracle instance of ChEMBL shall be used) e.g. by running: `pip install cx-oracle`.


Set the following environment variables:

`CHEMBL_DB=[DB_TYPE]://[USERNAME]:[PASSWORD]@[HOST]:[PORT]/[DB_NAME]` (Note that this is not a valid address and needs to be replaced by your ChEMBL connection.)
Example: "CHEMBL_DB=mysql://user:password@host:4233/chembl_30" or "CHEMBL_DB=oracle://..."   

`CHEMBL_VERSION=30`

! The current pipeline works with ChEMBL version 26 and more recent versions (but not older ones) !



## Getting Started

Run the pipeline with the following command:

`run-ot-pipeline-v2 genes.csv`

where `genes.csv` is a file with one Ensembl Gene ID per line in the first column of the file. 
No header is required. 
Any other present columns will be ignored. 



## Bucket assignment details for all modality workflows

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
Bucket 6 (PhosphoSitePlus, mUbiSiDa, Kim et al. 2011) Targets with reported ubiquitination sites in PhosphoSitePlus, mUbiSiDa (2013), or [Kim et al. 2011](https://www.sciencedirect.com/science/article/pii/S1097276511006757)  
Bucket 7 ([Mathieson et al. 2018](https://www.nature.com/articles/s41467-018-03106-1) and [Rolfs et al. 2021](https://www.nature.com/articles/s41467-021-26842-3)) Targets with available half-life data  
Bucket 8 (ChEMBL) Targets with SM in ChEMBL with activity on target-based assay (pChEMBL ≥5) 

PROTAC_location_Bucket  
In addition to PROTAC tractability buckets, there is an additional "PROTAC_location_Bucket", which allows to assess
whether a target's location is suitable for the PROTAC approach. UniProt and GO locations are used. 
This bucket does not have a binary evaluation; Instead the target location is assessed with a score ranging from 1 to 7, depending on the favorability of the 
location for a PROTAC approach and the confidence level provided by the evidence source (UniProt/GO).  
Here, good location = [cytoplasm or nucleus], grey location = [membrane, including cell membrane and organelle membranes], and bad location = anywhere else.  
The PROTAC_location_Bucket score is then defined by:  
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
Bucket 3 (ChEMBL) Targets with drug (other than SM or AB) in clinical phase 1  



## Provided data besides bucket assignments

In the primary TSV output the first columns provide general target information. 
This includes Ensembl gene ID and gene symbol, UniProt ID, Entry name and Protein names, 
information from the Illuminating the Druggable Genome (IDG) project with 'IDG_family' (family adherence), 
'IDG_dtoclass' (Drug Target Ontology class), 'IDG_tdl' (Target Development Level), and 'IDG_fam' (global IDG family classification), 
as well as GO Biological Process and Molecular Function terms.

There are also several columns following the bucket columns providing additional data, which can be exploited by the user.  

All workflows (*W*) have two bucket summary columns:  
'Bucket_sum_*W*' = sum of assessments for all binary buckets  
'Top_bucket_*W*' = highest scoring bucket (lowest bucket number with a positive assignment)  

Additionally, for all workflows overall tractability categories are assigned. 
The attribution of a target to certain bucket results in a workflow category assignment (see details below for each workflow), 
indicated with a score in a category column. 
The categories are also summarized per workflow in the 'Top_Category_*W*' column, containing the name of the target's highest assigned category or 'Unknown'.


#### Small molecule (SM) workflow
SM category assignment rules as 'Top_Category_sm' (in descending order):
- 'Clinical_Precedence_sm' is assigned if any of the SM clinical buckets 1-3 are positive. 
(The category score is calculated by 1 * ['Bucket_1_sm'] + 0.7 * ['Bucket_2_sm'] + 0.2 * ['Bucket_3_sm'].)
- 'Discovery_Precedence_sm' is assigned if SM bucket 4 or 5 is the highest assigned bucket. 
(The category score is calculated by 0.7 * ['Bucket_4_sm'] + 0.3 * ['Bucket_5_sm'].)
- 'Predicted_Tractable_sm' is assigned if SM bucket 6, 7 or 8 is the highest assigned bucket. 
(The category score is calculated by 0.7 * ['Bucket_6_sm'] + 0.3 * ['Bucket_7_sm'] + 0.3 * ['Bucket_8_sm'].)

Columns 'drug_chembl_ids_sm', 'drug_names_sm', and 'clinical_phases_sm' contain data related to clinical buckets 1-3.  
Column 'PDB_Known_Ligand' contains the PDB codes of structures referring to bucket 4 annotation, and column 'PDB_other' contains 
additional PDB codes where structures are available, but the ligand criteria for bucket 4 assignment are not met.  
Column 'High_Quality_ChEMBL_compounds' contains the amount of compounds found in ChEMBL related to bucket 5 assignment.  
Column 'DrugEBIlity_score' containes the score provided by the tool DrugEBIlity (related to assignment of buckets 6 and 7).   
Column 'Small_Molecule_Druggable_Genome_Member' contains 'Y' for Yes or 'N' for No related to bucket 8 assignment.  


#### Antibody (AB) workflow
AB category assignment rules as 'Top_Category_ab' (in descending order):
- 'Clinical_Precedence_ab' is assigned if any of the AB clinical buckets 1-3 are positive. 
(The category score is calculated by 1 * ['Bucket_1_ab'] + 0.7 * ['Bucket_2_ab'] + 0.2 * ['Bucket_3_ab'].)
- 'Predicted_Tractable_ab_High_confidence' is assigned if AB bucket 4 or 5 is the highest assigned bucket. 
(The category score is calculated by 0.7 * ['Bucket_4_ab'] + 0.3 * ['Bucket_5_ab'].)
- 'Predicted_Tractable_ab_Medium_to_low_confidence' is assigned if AB bucket 6, 7, 8 or 9 is the highest assigned bucket. 
(The category score is calculated by 0.4 * ['Bucket_6_ab'] + 0.25 * ['Bucket_7_ab'] + 0.25 * ['Bucket_8_ab'] + 0.1 * ['Bucket_9_ab'].)

Columns 'drug_chembl_ids_ab', 'drug_names_ab', and 'clinical_phases_ab' contain data related to clinical buckets 1-3.  
Columns 'Uniprot_high_conf_loc', 'GO_high_conf_loc', 'Uniprot_med_conf_loc', and 'GO_med_conf_loc' contain the subcellular 
location names and their respective evidence codes provided by the evidence source (UniProt/GO), data related to buckets 4, 5, 6 and 8.  
Confidence for location evidence codes (provided by UniProt):  
high_conf_evidence = 'ECO:0000269' or 'ECO:0000305'  
medium_conf_evidence = any other type of UniProt evidence  
(Details on UniProt evidence codes can be found here: https://www.uniprot.org/help/evidences)  
Confidence for location evidence types (provided by GO):  
evidence_types = {'EXP': 'High', 'IDA': 'High', 'IPI': 'High', 'TAS': 'High', 'IMP': 'High', 'IGI': 'High', 'IEP': 'High',
'ISS': 'Medium', 'ISO': 'Medium', 'ISA': 'Medium', 'ISM': 'Medium', 'IGC': 'Medium', 'IBA': 'Medium', 'IBD': 'Medium', 'IKR': 'Medium', 'IRD': 'Medium', 'RCA': 'Medium', 'IEA': 'Medium',
'NAS': 'Low', 'IC': 'Low', 'ND': 'Low', 'NR': 'Low'}  
(Details on GO evidence codes can be found here: http://geneontology.org/docs/guide-go-evidence-codes/)  
Columns 'Transmembrane' and 'Signal_peptide' contain the labels ['TRANSMEM'] (indicative for a transmembrane protein) and ['SIGNAL'] 
(indicative for a signal peptide), provided by UniProt and related to bucket 7 assignment.  
Column 'HPA_main_location' contains the main location data provided by HPA if it is labelled with "Approved", "Enhanced" or "Supported" reliability; data related to bucket 9. 
(The complete HPA data file can be found here: https://www.proteinatlas.org/download/subcellular_location.tsv.zip and 
more on HPA Reliability scores can be found here: https://www.proteinatlas.org/about/assays+annotation#if_reliability_score.)  


#### PROTAC workflow
PROTAC category assignment rules as 'Top_Category_PROTAC' (in descending order):
- 'Clinical_Precedence_PROTAC' is assigned if any of the clinical buckets 1-3 are positive.
- 'Literature_Precedence_PROTAC' is assigned if literature bucket 4 is positive.
- 'Discovery_Opportunity_PROTAC' is assigned for a given target if bucket (5 or 6) and bucket 7 and bucket 8 are assigned and 
the PROTAC_location_Bucket annotates high or medium confidence for "good" or "grey" location (location bucket score of 1, 2, 3 or 4).

Columns 'mentioned_in_PROTAC_literature', 'literature_count_PROTAC', 'pub_id', 'full_id' and 'title' relate to the fully automatic PROTAC literature 
assessment based on abstracts from EuropePMC, indicate whether a target has been identified in PROTAC literature and provide the related evidence. 
This data is not represented by a bucket.  
Column 'Uniprot_keyword' indicates whether a target has the keyword: "Ubl conjugation [KW-0832]" annotated in UniProt and relates to bucket 5
Column 'Uniprot_PTM' contains additional information on posttranslational modifications from UniProt and is not associated with a bucket.  
Column 'Uniprot_CrossLink' contains additional information on crosslinks from UniProt and is not associated with a bucket.  
Columns 'Ub_PhosphoSitePlus', 'Ub_mUbiSiDa_2013' and 'number_of_ubiquitination_sites' relate to bucket 6 and represent information available 
from the two databases PhosphoSitePlus and mUbiSiDa and the dataset from Kim et al., respectively.  
Columns 'Max_halflife' and 'Min_halflife' contain the max/min values of half-life data means measured in different cell types from Mathieson et al. and Rolfs et al. and provides more detail to bucket 7 assessment. 
The available data is provided in columns 'Bcell_mean', 'NKcell_mean', 'Hepatocytes_mean', 'MouseNeuorons_mean'.  
Column 'count_compound_chembl_ids_PROTAC' contains the number of compounds found in ChEMBL related to bucket 8.

##### Additional files & tools:  
Files from the automated assessment process where specific targets are first tagged using an implementation from Europe PMC (Europe PMC Annotations API provides text mining annotations contained in abstracts and open access full text articles, 
using the W3C Open Annotation Data Model https://europepmc.org/AnnotationsApi), where we consider annotations from abstracts, followed by mapping to their respective human UniProt IDs. This enables detailed target information to be provided to 
the user in the form of three different lists (that can be found in the pipelines "fetched data" output directory): 1) a complete list with all detected targets (named "protac_literature_evidence.csv") and two subsets 2) targets that are 
associated with the ubiquitin ligase complex, including E3 ligases (filtering is based on evidence from the respective UniProt entry text; named "protac_literature_evidence_ubi_ligase_complex.csv"), and 3) the remaining targets (named 
"protac_literature_evidence_filtered.csv"). A literature count is provided per target with publication titles, links and full abstract texts.  

To aid further in literature data processing with respect to detection of PROTAC relevant information additional Natural Language Processing (NLP) is performed with a Named Entity Recognition (NER) model that has been trained and implemented for 
automatic detection of the entities "PROTAC_TARGET", "PROTAC_NAME", and "E3_LIGASE". The ScispaCy ([Neumann et al. 2019](https://aclanthology.org/W19-5034/)) pre-trained model "en_ner_jnlpba_md" (trained on the JNLPBA corpus) was used to retrain the PROTAC-NER model on 164 abstracts 
that mention PROTACs, with 427 manually labelled entities (280 "PROTAC_TARGET", 57 "PROTAC_NAME", and 90 "E3_LIGASE"). The PROTAC-NER abstract screening is performed on all PROTAC literature abstracts when running the pipeline and the results are 
provided to the user as addition to the literature evidence list. The trained PROTAC-NER model is also provided with the pipeline code and it can be used with spaCy to detect PROTAC entities in any text.  


#### Other clinical evidence workflow
OtherCl category assignment rule:
- 'Clinical_Precedence_othercl' is assigned if any of the OtherCl clinical buckets 1-3 are positive.
