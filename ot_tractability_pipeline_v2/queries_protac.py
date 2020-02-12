# SQLqueries for PROTACs
from ot_tractability_pipeline_v2.queries_sm import CHEMBL_VERSION

chembl_clinical_protac = '''
SELECT DISTINCT mh.parent_molregno, di.efo_id, di.efo_term, di.max_phase_for_ind
FROM {0}.molecule_dictionary md,
	{0}.molecule_hierarchy mh,
	{0}.drug_indication di
WHERE md.molregno = mh.molregno
AND md.molregno = di.molregno

AND md.molecule_type = 'PROTAC'
'''.format(CHEMBL_VERSION)

chembl_clinical_protac_targets = '''
SELECT DISTINCT mh.parent_molregno, 
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
'''.format(CHEMBL_VERSION)