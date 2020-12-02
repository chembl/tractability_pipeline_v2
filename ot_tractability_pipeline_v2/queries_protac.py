# SQLqueries for PROTACs
from ot_tractability_pipeline_v2.queries_sm import CHEMBL_VERSION

chembl_clinical_protac = '''
    SELECT DISTINCT mh.parent_molregno, 
        di.efo_id, 
        di.efo_term, 
        di.max_phase_for_ind
    FROM {0}.molecule_dictionary md,
    	{0}.molecule_hierarchy mh,
    	{0}.drug_indication di
    WHERE md.molregno = mh.molregno
    AND md.molregno = di.molregno
    
    AND md.molecule_type = 'PROTAC'
    '''.format(CHEMBL_VERSION)

chembl_clinical_protac_targets = '''
    SELECT DISTINCT mh.parent_molregno, 
        md.chembl_id AS drug_chembl_id,
    	md.pref_name AS drug_name, 
    	md.max_phase,
    	td.tid,
    	td.pref_name AS target_name, 
    	td.target_type,
    	cs.accession AS accession,
    	dm.action_type AS moa_chembl,
    	mr.ref_type,
    	mr.ref_id,
    	mr.ref_url,
    	dm.site_id
    FROM {0}.molecule_dictionary md,
    	{0}.molecule_hierarchy mh,
    	{0}.drug_mechanism dm,
    	{0}.target_dictionary td,
    	{0}.target_components tc,
    	{0}.component_sequences cs,
    	{0}.mechanism_refs mr
    WHERE md.molregno = dm.molregno
    AND md.molregno = mh.molregno
    AND dm.tid = td.tid
    AND td.tid = tc.tid
    AND tc.component_id = cs.component_id
    AND dm.mec_id = mr.mec_id
    
    AND md.molecule_type = 'PROTAC'
    AND td.tax_id = 9606
    AND td.target_type LIKE '%PROTEIN%'
    '''.format(CHEMBL_VERSION)

# =============================================================================
# get all small molecules from ChEMBL
# =============================================================================
chembl_small_mol = """
    SELECT DISTINCT mh.parent_molregno, 
         md.molregno, 
         md.pref_name AS compound_name
    FROM {0}.molecule_dictionary md,
        {0}.molecule_hierarchy mh
    WHERE md.molregno = mh.molregno
    
    AND md.molecule_type = 'Small molecule'
    """.format(CHEMBL_VERSION)


# =============================================================================
# get all small molecules and their targets with a pChEMBL >=5
# =============================================================================
chembl_small_mol_active_targets = """
    SELECT DISTINCT mh.parent_molregno,
        md.chembl_id AS compound_chembl_id,
        md.pref_name AS compound_name,
        act.pchembl_value,
        td.tid,
        td.chembl_id AS target_chembl_id,
        td.pref_name AS target_name,
        td.target_type,
        cs.accession AS accession
    FROM {0}.molecule_dictionary md,
        {0}.molecule_hierarchy mh,
        {0}.activities act,
        {0}.assays a,
        {0}.target_dictionary td,
        {0}.target_components tc,
        {0}.component_sequences cs,
        {0}.compound_properties cp
    WHERE md.molregno = mh.molregno
    AND md.molregno = act.molregno 
    AND md.molregno = cp.molregno
    AND act.assay_id = a.assay_id
    AND a.tid = td.tid
    AND td.tid = tc.tid
    AND tc.component_id = cs.component_id
    AND td.tax_id = 9606
    AND td.target_type LIKE '%PROTEIN%'
    AND act.pchembl_value >= 5
    AND md.molecule_type not in ('Oligonucleotide', 'Oligosaccharide')
    AND ((md.molecule_type in ('Protein', 'Unclassified', 'Unknown') AND cp.mw_freebase <= 1500) OR (md.molecule_type = 'Small molecule'))
    """.format(CHEMBL_VERSION)
#    AND md.molecule_type = 'Small molecule'

    