from __future__ import print_function, absolute_import, division

from six import iteritems
from warnings import warn
import json

import pandas as pd

import cobrame


# MetabolicReactions present in iOL1650 that have been removed for iLE1678
removed_reactions = ['ALPATG160pp1', 'ALPATE160pp1', 'ATPM',
                     'CITLY-CPLX_2tpr3dpcoa', 'PFL_act']

# PFL enzymes are activated by a glycl radical group
pfl_isozymes = ['PYRUVFORMLY-INACTIVE-CPLX',
                'PYRUVFORMLY-MONOMER_EG11784-MONOMER',
                'KETOBUTFORMLY-INACT-MONOMER',
                'EG11910-MONOMER_dimer_EG11911-MONOMER']


def update_metabolite_formulas(m_model):
    # Formulas of metabolites not included in metabolites.txt
    # 3a1hac1p_c has a typo in the formula in metabolites.txt
    formulas = [('4fe4s_c', 'Fe4S4'), ('2fe2s_c', 'Fe2S2'),
                ('LI_c', 'Li'), ('3a1hac1p_c', 'C3H7NO5P'),
                ]

    for met, formula in formulas:
        try:
            met_obj = m_model.metabolites.get_by_id(met)
        except KeyError:
            warn('Creating new metabolite (%s)' % met)
            met_obj = cobrame.Metabolite(met)
            m_model.add_metabolites([met_obj])
        met_obj.formula = formula


def correct_complex_modifications(model):
    # Known modification (w/ high confidence)
    mod = model.process_data.mod_acetyl_c
    mod.stoichiometry = {'accoa_c': -1, 'coa_c': 1}
    mod._element_contribution = {'C': 2, 'H': 2, 'O': 1}

    # Per 23597401, this is a possible stoichiometry for this modification
    mod = model.process_data.mod_NiFeCoCN2_c
    mod.stoichiometry = {'fe2_c': -1, 'co2_c': -1, 'nadh_c': -2, 'nad_c': 2,
                         'cbp_c': -2, 'atp_c': -4, 'amp_c': 2, 'adp_c': 2,
                         'ppi_c': 2, 'pi_c': 4, 'h_c': 2, 'h2o_c': -1,
                         'ni2_c': -1}
    mod._element_contribution = {'C': 3, 'Fe': 1, 'N': 2, 'Ni': 1, 'O': 1}

    # Updates to mod_lipo stoichiometry
    model.process_data.mod_lipo_c.stoichiometry['h_c'] = 2
    model.process_data.mod_lipo_c_alt.stoichiometry['h_c'] = 1

    # PFL is activated by a glycl group by PFLACTENZ-MONOMER
    mod = model.process_data.mod_glycl_c
    mod.enzyme = ['PFLACTENZ-MONOMER', 'FLAVODOXIN1-MONOMER']
    mod.stoichiometry = {'FLAVODOXIN1-MONOMER': -1, 'amet_c': -1,
                         'dad__2_c': 1, 'FLAVODOXIN1-MONOMER_mod_Oxidized': 1,
                         'met__L_c': 1}
    mod._element_contribution = {'H': -1}


def correct_reaction_matrix(reaction_matrix_dict):
    for r in removed_reactions:
        reaction_matrix_dict.pop(r)

    # Per 10510271, reaction to reduce CU(II) to CU(I)
    reaction_matrix_dict['CU2R'] = {'cu2_c': -1, 'cu_c': 1, 'nadh_c': -1,
                                    'nad_c': 1, 'h_c': 1}
    return reaction_matrix_dict


def correct_reaction_info_frame(df):
    for r in removed_reactions:
        df = df.drop(r)

    # These reactions are irreversible and are highly active if not corrected
    df.loc['PPKr', 'is_reversible'] = 0
    df.loc['PPK2r', 'is_reversible'] = 0

    # Per 10510271, reaction to reduce CU(II) to CU(I)
    df = df.append(pd.Series({'description': 'CU2 Reduction',
                              'is_reversible': 0, 'is_spontaneous': 0},
                             name='CU2R'))

    return df


def correct_reaction_stoichiometries(model, file_name):
    df = pd.read_excel(file_name, index_col=0)
    for d in df.index:
        stoich_data = model.process_data.get_by_id(d)
        stoich_dict = json.loads(df.loc[d, 'Stoich Change'].replace("'", "\""))
        for met, stoich in iteritems(stoich_dict):
            stoich_data.stoichiometry[met] = stoich

        stoich_data._update_parent_reactions()


def correct_enzyme_reaction_association_frame(df):

    # This ATPS4rpp should be associated with the modified version of this
    # complex
    df = \
        df.applymap(lambda x: x.replace('ATPSYN-CPLX_EG10106-MONOMER',
                                        'ATPSYN-CPLX_EG10106-MONOMER_mod_mg2'))

    df = df.applymap(
        lambda x: x.replace('EG11910-MONOMER_dimer',
                            'EG11910-MONOMER_dimer_EG11911-MONOMER'))
    for cplx in pfl_isozymes:
        df = df.applymap(
            lambda x: x.replace(cplx, cplx + '_mod_glycl'))

    # Per 10510271, reaction to reduce CU(II) to CU(I)
    df.loc['CU2R', 'Complexes'] = 'NADH-DHII-MONOMER_mod_mg2_mod_cu_mod_fad'

    return df


def correct_trna_modifications(mod):
    """
    Apply corrections to tRNA modification procedures
    """
    # Per: 23543739, grxD is involved in repairing miaB 4fe4s prosthetic group
    # after accepting electron
    correct_mod = mod['ms2i6A_at_37']['carriers']
    correct_mod["CPLX0-7817"] = correct_mod.pop('fldrd_c')
    correct_mod["CPLX0-7817_mod_Oxidized"] = correct_mod.pop('fldox_c')

    # per 24914049, the 4fe4s prosthetic group in YdaO helps donate sulfur,
    # an additional ferredoxin is not required
    correct_mod = mod['s2C_at_32']['carriers']
    correct_mod.pop('trdrd_c')
    correct_mod.pop('trdox_c')
    correct_mod = mod['s2C_at_32']['metabolites']
    correct_mod['h_c'] = -1
    mod['s2C_at_32']['machines'] = ['YdaO_dim_mod_4fe4s']

    # Per PMID 10753862, mechanism does not require thioredoxin or other
    # reducing equivalent
    correct_mod = mod['s4U_at_8']['carriers']
    correct_mod.pop('trdrd_c')
    correct_mod.pop('trdox_c')
    correct_mod = mod['s4U_at_8']['metabolites']
    correct_mod['h_c'] = -1

    # Assume thioredoxin 1 though others could perform this role
    correct_mod = mod['mnm5s2U_at_34']['carriers']
    correct_mod['RED-THIOREDOXIN-MONOMER'] = correct_mod.pop('trdrd_c')
    correct_mod['RED-THIOREDOXIN-MONOMER_mod_Oxidized'] = correct_mod.pop(
        'trdox_c')

    # Per 19322199, atp is required for this modification
    correct_mod = mod['ac4C_at_34']['metabolites']
    correct_mod['atp_c'] = -1
    correct_mod['h2o_c'] = -1
    correct_mod['h_c'] = 1
    correct_mod['adp_c'] = 1
    correct_mod['pi_c'] = 1

    correct_mod = mod['m6t6A_at_37']['metabolites']
    correct_mod['h_c'] = 1

    # Cobalamin stimulates activity of QueG, but is not required,
    # per PMID:21502530
    mod['Q_at_34']['machines'] = ['Tgt_hexa_mod_6:zn2', 'QueA_mono',
                                  'QueG_mono']
    correct_mod = mod['Q_at_34']['metabolites']
    correct_mod['h2o_c'] = 1
    correct_mod['h_c'] = 1
    mod['Q_at_34']['carriers'] = {}
    correct_mod = mod['Q_at_34']['carriers']
    correct_mod['RED-THIOREDOXIN-MONOMER'] = - 1
    correct_mod['RED-THIOREDOXIN-MONOMER_mod_Oxidized'] = 1

    # Changed stoichiometry per PMIDs: 23676670, 25855808, 26681692.
    # Reaction forming 5-hydroxyuridine is not known
    mod['cmo5U_at_34']['metabolites'] = {'phpyr_c': 1, 'pphn_c': -1,
                                         'amet_c': -2, 'h2o_c': 1, 'h_c': 1,
                                         'ahcys_c': 2}
    return mod


def correct_rrna_modifications(mod):
    mod['m7G_at_2069']['machine'] = 'RlmL_dim'
    return mod


def correct_complex_stoichiometry(stoichiometry):
    # add in missing entries
    stoichiometry["CPLX0-7617"] = {"protein_b0156": 2}

    # Complex Updated based on iJO1366 GPR
    stoichiometry['EG11910-MONOMER_dimer_EG11911-MONOMER'] = \
        {'EG11910-MONOMER_dimer': 1, 'EG11911-MONOMER': 1}

    # should actually be a dimer PMID 24914049
    stoichiometry.pop('YdaO_mono')
    stoichiometry['YdaO_dim'] = {"protein_b1344": 2}

    # Replaced by lipoprotein reactions from PMID:25227965
    stoichiometry.pop('EG10544-MONOMER')

    stoichiometry['ATPSYN-CPLX_EG10106-MONOMER'] = {'protein_b3731': 1.0,
                                                    'protein_b3732': 3.0,
                                                    'protein_b3733': 1.0,
                                                    'protein_b3734': 3.0,
                                                    'protein_b3735': 1.0,
                                                    'protein_b3736': 2.0,
                                                    'protein_b3737': 10.0,
                                                    'protein_b3738': 1.0,
                                                    'protein_b3739': 1.0}

    return stoichiometry


def correct_complex_modification_dict(stoichiometry):
    # specific patches. Not used in iOL1650 but included as a complex in
    # protein_modification.txt
    for i in ['CPLX0-246_CPLX0-1342_mod_1:SH',
              'EG10544-MONOMER_mod_palmitate', 'EG11597-MONOMER_mod_coo']:
        stoichiometry.pop(i)

    stoichiometry['CPLX0-246_CPLX0-1342_mod_pydx5p'] = {'core_enzyme':
                                                        'CPLX0-246_CPLX0-1342',
                                                        'modifications': {
                                                           'pydx5p_c': 1}}
    stoichiometry["IscS_mod_2:pydx5p"] = {"core_enzyme": "IscS",
                                          "modifications": {"pydx5p_c": 2}}
    stoichiometry["YdaO_dim_mod_4fe4s"] = {"core_enzyme": "YdaO_dim",
                                           "modifications": {"4fe4s_c": 1}}

    # new complex for ATP synthase complex with atpI subunit
    stoichiometry['ATPSYN-CPLX_EG10106-MONOMER_mod_mg2'] = \
        {'core_enzyme': 'ATPSYN-CPLX_EG10106-MONOMER',
         'modifications': {'mg2_c': 1.}}

    for cplx in pfl_isozymes:
        stoichiometry[cplx + '_mod_glycl'] = {'core_enzyme': cplx,
                                              "modifications": {'glycl_c': 1}}

    return stoichiometry
