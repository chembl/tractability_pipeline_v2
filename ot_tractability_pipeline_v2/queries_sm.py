# SQL queries for small molecules
# The length of the SQL queries makes the Python code far too long, moved to separate file

import os

try:
    chembl_number = int(os.getenv('CHEMBL_VERSION'))
    CHEMBL_VERSION = "CHEMBL_{}".format(chembl_number)
except TypeError:
    print('''
    Please set the CHEMBL version '''
          )
    raise

            
chembl_clinical_small_mol = """
    SELECT DISTINCT mh.parent_molregno, 
        md.molregno, 
        md.pref_name,
        di.efo_id, 
        di.efo_term, 
        di.max_phase_for_ind
    FROM {0}.molecule_dictionary md,
    	{0}.molecule_hierarchy mh,
    	{0}.drug_indication di
    WHERE md.molregno = mh.molregno
    AND md.molregno = di.molregno
    
    AND md.molecule_type = 'Small molecule'
        """.format(CHEMBL_VERSION)

chembl_clinical_targets = """
    SELECT DISTINCT mh.parent_molregno, 
        md.chembl_id AS drug_chembl_id,
        md.pref_name AS drug_name, 
        md.max_phase,
        td.tid,
        td.chembl_id AS target_chembl_id,
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
    
    AND md.molecule_type = 'Small molecule'
    AND td.tax_id = 9606
    AND td.target_type LIKE '%PROTEIN%'
    """.format(CHEMBL_VERSION)

pchembl_q = '''
    SELECT DISTINCT td.chembl_id AS target_chembl_id, 
        td.tid,
        cs.accession, 
        pf.protein_class_desc, 
        md.molregno, 
        md.chembl_id AS compound_chembl_id,
        cst.canonical_smiles,
        cp.cx_logd,
        cp.aromatic_rings,
        cp.qed_weighted,
        d.year
    FROM {0}.target_dictionary td,
        {0}.target_components tc,
        {0}.component_sequences cs,
        {0}.component_class cc,
        {0}.protein_classification pf,
        {0}.assays ass,
        {0}.activities act,
        {0}.molecule_dictionary md,
        {0}.compound_structures cst,
        {0}.compound_properties cp,
        {0}.docs d
    WHERE td.tid = tc.tid
    AND tc.component_id = cs.component_id
    AND cs.component_id = cc.component_id
    AND cc.protein_class_id = pf.protein_class_id
    AND td.tid = ass.tid
    AND ass.assay_id = act.assay_id
    AND act.molregno = md.molregno
    AND md.molregno = cst.molregno
    AND md.molregno = cp.molregno
    AND act.doc_id = d.doc_id
    AND td.tax_id = 9606
    AND td.target_type LIKE '%PROTEIN%'
    AND ass.assay_type = 'B'
    AND ass.relationship_type = 'D'
    AND (act.data_validity_comment is NULL OR act.data_validity_comment = 'Manually validated')
    AND act.potential_duplicate = 0
    AND (act.pchembl_value is not NULL AND pchembl_value >= 5.5)
    AND md.molecule_type not in ('Oligonucleotide', 'Oligosaccharide')
    AND ((md.molecule_type in ('Protein', 'Unclassified', 'Unknown') AND cp.mw_freebase <= 1500) OR (md.molecule_type = 'Small molecule'))
    '''.format(CHEMBL_VERSION)

nm_q = '''
    SELECT DISTINCT td.chembl_id AS target_chembl_id, 
        td.tid, 
        cs.accession, 
        pf.protein_class_desc, 
        md.molregno, 
        md.chembl_id AS compound_chembl_id,
        cst.canonical_smiles,
        cp.cx_logd,
        cp.aromatic_rings,
        cp.qed_weighted,
        d.year
    FROM {0}.target_dictionary td,
        {0}.target_components tc,
        {0}.component_sequences cs,
        {0}.component_class cc,
        {0}.protein_classification pf,
        {0}.assays ass,
        {0}.activities act,
        {0}.molecule_dictionary md,
        {0}.compound_structures cst,
        {0}.compound_properties cp,
        {0}.docs d
    WHERE td.tid = tc.tid
    AND tc.component_id = cs.component_id
    AND cs.component_id = cc.component_id
    AND cc.protein_class_id = pf.protein_class_id
    AND td.tid = ass.tid
    AND ass.assay_id = act.assay_id
    AND act.molregno = md.molregno
    AND md.molregno = cst.molregno
    AND md.molregno = cp.molregno
    AND act.doc_id = d.doc_id
    AND td.tax_id = 9606
    AND td.target_type LIKE '%PROTEIN%'
    AND ass.assay_type = 'B'
    AND ass.relationship_type = 'D'
    AND (act.data_validity_comment is NULL OR act.data_validity_comment = 'Manually validated')
    AND act.potential_duplicate = 0
    AND (act.pchembl_value is NULL AND act.standard_type in ('AC50', 'CC50', 'EC50', 'GI50', 'IC50', 'IC90', 'IC95', 'IC99', 'Kd', 'Ki', 'LC50', 'MIC', 'MIC50', 'Potency', 'Kinact', 'KB', 'Activity', 'Ke', 'KA', 'IC100') AND act.standard_units = 'nM' AND act.standard_value <= 3000 AND act.standard_relation in ('<', '='))
    AND md.molecule_type not in ('Oligonucleotide', 'Oligosaccharide')
    AND ((md.molecule_type in ('Protein', 'Unclassified', 'Unknown') AND cp.mw_freebase <= 1500) OR (md.molecule_type = 'Small molecule'))
    '''.format(CHEMBL_VERSION)

km_kon_q = '''
    SELECT DISTINCT td.chembl_id AS target_chembl_id, 
        td.tid, 
        cs.accession, 
        pf.protein_class_desc, 
        md.molregno, 
        md.chembl_id AS compound_chembl_id,
        cst.canonical_smiles,
        cp.cx_logd,
        cp.aromatic_rings,
        cp.qed_weighted,
        d.year
    FROM {0}.target_dictionary td,
        {0}.target_components tc,
        {0}.component_sequences cs,
        {0}.component_class cc,
        {0}.protein_classification pf,
        {0}.assays ass,
        {0}.activities act,
        {0}.molecule_dictionary md,
        {0}.compound_structures cst,
        {0}.compound_properties cp,
        {0}.docs d
    WHERE td.tid = tc.tid
    AND tc.component_id = cs.component_id
    AND cs.component_id = cc.component_id
    AND cc.protein_class_id = pf.protein_class_id
    AND td.tid = ass.tid
    AND ass.assay_id = act.assay_id
    AND act.molregno = md.molregno
    AND md.molregno = cst.molregno
    AND md.molregno = cp.molregno
    AND act.doc_id = d.doc_id
    AND td.tax_id = 9606
    AND td.target_type LIKE '%PROTEIN%'
    AND ass.assay_type = 'B'
    AND ass.relationship_type = 'D'
    AND (act.data_validity_comment is NULL OR act.data_validity_comment = 'Manually validated')
    AND act.potential_duplicate = 0
    AND (act.pchembl_value is NULL AND act.standard_type in ('k_on', 'Km'))
    AND md.molecule_type not in ('Oligonucleotide', 'Oligosaccharide')
    AND ((md.molecule_type in ('Protein', 'Unclassified', 'Unknown') AND cp.mw_freebase <= 1500) OR (md.molecule_type = 'Small molecule'))
    '''.format(CHEMBL_VERSION)

D_Tm_q = '''
    SELECT DISTINCT td.chembl_id AS target_chembl_id, 
        td.tid, 
        cs.accession, 
        pf.protein_class_desc, 
        md.molregno, 
        md.chembl_id AS compound_chembl_id,
        cst.canonical_smiles,
        cp.cx_logd,
        cp.aromatic_rings,
        cp.qed_weighted,
        d.year
    FROM {0}.target_dictionary td,
        {0}.target_components tc,
        {0}.component_sequences cs,
        {0}.component_class cc,
        {0}.protein_classification pf,
        {0}.assays ass,
        {0}.activities act,
        {0}.molecule_dictionary md,
        {0}.compound_structures cst,
        {0}.compound_properties cp,
        {0}.docs d
    WHERE td.tid = tc.tid
    AND tc.component_id = cs.component_id
    AND cs.component_id = cc.component_id
    AND cc.protein_class_id = pf.protein_class_id
    AND td.tid = ass.tid
    AND ass.assay_id = act.assay_id
    AND act.molregno = md.molregno
    AND md.molregno = cst.molregno
    AND md.molregno = cp.molregno
    AND act.doc_id = d.doc_id
    AND td.tax_id = 9606
    AND td.target_type LIKE '%PROTEIN%'
    AND ass.assay_type = 'B'
    AND ass.relationship_type = 'D'
    AND (act.data_validity_comment is NULL OR act.data_validity_comment = 'Manually validated')
    AND act.potential_duplicate = 0
    AND (act.pchembl_value is NULL AND act.standard_type in ('deltaTm', 'Tm') AND act.standard_units = 'degrees C' AND act.standard_value >= 2 AND act.standard_relation in ('>', '='))
    AND md.molecule_type not in ('Oligonucleotide', 'Oligosaccharide')
    AND ((md.molecule_type in ('Protein', 'Unclassified', 'Unknown') AND cp.mw_freebase <= 1500) OR (md.molecule_type = 'Small molecule'))
    '''.format(CHEMBL_VERSION)

residual_act_q = '''
    SELECT DISTINCT td.chembl_id AS target_chembl_id, 
        td.tid, 
        cs.accession, 
        pf.protein_class_desc, 
        md.molregno, 
        md.chembl_id AS compound_chembl_id,
        cst.canonical_smiles,
        cp.cx_logd,
        cp.aromatic_rings,
        cp.qed_weighted,
        d.year
    FROM {0}.target_dictionary td,
        {0}.target_components tc,
        {0}.component_sequences cs,
        {0}.component_class cc,
        {0}.protein_classification pf,
        {0}.assays ass,
        {0}.activities act,
        {0}.molecule_dictionary md,
        {0}.compound_structures cst,
        {0}.compound_properties cp,
        {0}.docs d
    WHERE td.tid = tc.tid
    AND tc.component_id = cs.component_id
    AND cs.component_id = cc.component_id
    AND cc.protein_class_id = pf.protein_class_id
    AND td.tid = ass.tid
    AND ass.assay_id = act.assay_id
    AND act.molregno = md.molregno
    AND md.molregno = cst.molregno
    AND md.molregno = cp.molregno
    AND act.doc_id = d.doc_id
    AND td.tax_id = 9606
    AND td.target_type LIKE '%PROTEIN%'
    AND ass.assay_type = 'B'
    AND ass.relationship_type = 'D'
    AND (act.data_validity_comment is NULL OR act.data_validity_comment = 'Manually validated')
    AND act.potential_duplicate = 0
    AND (act.pchembl_value is NULL AND act.standard_type = 'Residual activity' AND act.standard_units = '%' AND act.standard_value <= 10 AND act.standard_relation in ('<', '='))
    AND md.molecule_type not in ('Oligonucleotide', 'Oligosaccharide')
    AND ((md.molecule_type in ('Protein', 'Unclassified', 'Unknown') AND cp.mw_freebase <= 1500) OR (md.molecule_type = 'Small molecule'))
    '''.format(CHEMBL_VERSION)

Imax_q = '''
    SELECT DISTINCT td.chembl_id AS target_chembl_id, 
        td.tid, 
        cs.accession, 
        pf.protein_class_desc, 
        md.molregno, 
        md.chembl_id AS compound_chembl_id,
        cst.canonical_smiles,
        cp.cx_logd,
        cp.aromatic_rings,
        cp.qed_weighted,
        d.year
    FROM {0}.target_dictionary td,
        {0}.target_components tc,
        {0}.component_sequences cs,
        {0}.component_class cc,
        {0}.protein_classification pf,
        {0}.assays ass,
        {0}.activities act,
        {0}.molecule_dictionary md,
        {0}.compound_structures cst,
        {0}.compound_properties cp,
        {0}.docs d
    WHERE td.tid = tc.tid
    AND tc.component_id = cs.component_id
    AND cs.component_id = cc.component_id
    AND cc.protein_class_id = pf.protein_class_id
    AND td.tid = ass.tid
    AND ass.assay_id = act.assay_id
    AND act.molregno = md.molregno
    AND md.molregno = cst.molregno
    AND md.molregno = cp.molregno
    AND act.doc_id = d.doc_id
    AND td.tax_id = 9606
    AND td.target_type LIKE '%PROTEIN%'
    AND ass.assay_type = 'B'
    AND ass.relationship_type = 'D'
    AND (act.data_validity_comment is NULL OR act.data_validity_comment = 'Manually validated')
    AND act.potential_duplicate = 0
    AND (act.pchembl_value is NULL AND act.standard_type in ('Activity', 'Imax') AND act.standard_units = '%' AND act.standard_value >= 70 AND act.standard_relation in ('>', '='))
    AND md.molecule_type not in ('Oligonucleotide', 'Oligosaccharide')
    AND ((md.molecule_type in ('Protein', 'Unclassified', 'Unknown') AND cp.mw_freebase <= 1500) OR (md.molecule_type = 'Small molecule'))
    '''.format(CHEMBL_VERSION)

Emax_q = '''
    SELECT DISTINCT td.chembl_id AS target_chembl_id, 
        td.tid, 
        cs.accession, 
        pf.protein_class_desc, 
        md.molregno, 
        md.chembl_id AS compound_chembl_id,
        cst.canonical_smiles,
        cp.cx_logd,
        cp.aromatic_rings,
        cp.qed_weighted,
        d.year
    FROM {0}.target_dictionary td,
        {0}.target_components tc,
        {0}.component_sequences cs,
        {0}.component_class cc,
        {0}.protein_classification pf,
        {0}.assays ass,
        {0}.activities act,
        {0}.molecule_dictionary md,
        {0}.compound_structures cst,
        {0}.compound_properties cp,
        {0}.docs d
    WHERE td.tid = tc.tid
    AND tc.component_id = cs.component_id
    AND cs.component_id = cc.component_id
    AND cc.protein_class_id = pf.protein_class_id
    AND td.tid = ass.tid
    AND ass.assay_id = act.assay_id
    AND act.molregno = md.molregno
    AND md.molregno = cst.molregno
    AND md.molregno = cp.molregno
    AND act.doc_id = d.doc_id
    AND td.tax_id = 9606
    AND td.target_type LIKE '%PROTEIN%'
    AND ass.assay_type = 'B'
    AND ass.relationship_type = 'D'
    AND (act.data_validity_comment is NULL OR act.data_validity_comment = 'Manually validated')
    AND act.potential_duplicate = 0
    AND (act.pchembl_value is NULL AND act.standard_type in ('Emax', 'Efficacy') AND act.standard_units = '%' AND act.standard_value >= 120 AND act.standard_relation in ('>', '='))
    AND md.molecule_type not in ('Oligonucleotide', 'Oligosaccharide')
    AND ((md.molecule_type in ('Protein', 'Unclassified', 'Unknown') AND cp.mw_freebase <= 1500) OR (md.molecule_type = 'Small molecule'))
    '''.format(CHEMBL_VERSION)
