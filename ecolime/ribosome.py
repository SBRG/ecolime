from cobrame import ComplexData, TranscribedGene, ModificationData
from cobrame.util.building import add_modification_data

import cobra

# Dictionary of {formation_step:[{metabolite:stoichiometry}]}
# Positive for reactants negative for products (complex formation convention)

# Leaving out Tig_mono trigger factor for now. It's a chaperone not part of formation

###### From Teddy's notes in original ME code #######
# What I'm trying to accomplish here: (old, but still useful)
# 1) 1 30Sp + 1 generic_16S -> 1 rib_30
# 2) 1 50Sp + 1 generic_23S + 1 generic_5S -> 1 rib_50
# 3) 1 rib_30 + 1 rib_50 -> 1 rib_70
# 4) 1 rib_70 + 1 b0884_assumedMonomer + 1 b1718_uniprotComplex  --> 1 rib_30_IF1_IF3 + 1 rib_50
# 5) 1 b3168_assumedMonomer_gtp + 1 rib_30_IF1_IF3 --> 1 rib_30_ini


# Mod:
#   Era_dim (assembly factor) + 2 gtp +2 h20->
#   30S_assembly_factor_gtp_hydrolying_assembly_phase_1_gtp -> 2 gtp + 2 h + 2 pi


# TODO Check how 2 gtp is added
ribosome_stoich = {'30_S_assembly_1_(215)': {'stoich': {'RpsD_mono': 1,
                                                        'RpsE_mono': 1,
                                                        'RpsF_mono': 1,
                                                        'RpsG_mono': 1,
                                                        'RpsH_mono': 1,
                                                        'RpsI_mono': 1,
                                                        'RpsK_mono': 1,
                                                        'RpsL_mono': 1,
                                                        'RpsM_mono': 1,
                                                        'RpsO_mono': 1,
                                                        'RpsP_mono': 1,
                                                        'RpsQ_mono': 1,
                                                        'RpsR_mono': 1,
                                                        'RpsS_mono': 1,
                                                        'RpsT_mono': 1,
                                                        'generic_16s_rRNAs': 1}},
                   '30_S_assembly_2_(21S)': {'stoich': {'mg2_c': 60,
                                                        'RpsA_mono': 1,
                                                        'RpsB_mono': 1,
                                                        'RpsC_mono': 1,
                                                        'RpsJ_mono': 1, # TODO is the RplJ stoich right?
                                                        'RpsN_mono': 1,
                                                        'RpsU_mono': 1,
                                                        'Sra_mono': 1}},
                   '50_S_assembly_1': {'stoich': {'generic_23s_rRNAs': 1,
                                                  'generic_5s_rRNAs': 1,
                                                  'RplA_mono': 1,
                                                  'RplB_mono': 1,
                                                  'RplC_mono': 1,
                                                  'RplD_mono': 1,
                                                  'RplE_mono': 1,
                                                  'RplI_mono': 1,
                                                  'RplJ_mono': 1,
                                                  'RplK_mono': 1,
                                                  'RplM_mono': 1,
                                                  'RplQ_mono': 1,
                                                  'RplS_mono': 1,
                                                  'RplT_mono': 1,
                                                  'RplU_mono': 1,
                                                  'RplV_mono': 1,
                                                  'RplW_mono': 1,
                                                  'RplX_mono': 1,
                                                  'RpmC_mono': 1,
                                                  'RpmG_mono': 1,
                                                  'RpmH_mono': 1,
                                                  'rpL7/12_mod_1:acetyl': 2}},
                   '50_S_assembly_2': {'stoich': {'mg2_c': 111,
                                                  'RplF_mono': 1,
                                                  'RplN_mono': 1,
                                                  'RplO_mono': 1,
                                                  'RplP_mono': 1,
                                                  'RplR_mono': 1,
                                                  'RplY_mono': 1,
                                                  'RpmA_mono': 1,
                                                  'RpmB_mono': 1,
                                                  'RpmD_mono': 1,
                                                  'RpmE_mono': 1,
                                                  'RpmF_mono': 1,
                                                  'RpmI_mono': 1,
                                                  'RpmJ_mono': 1,
                                                  'Tig_mono': 1}, # Leave Tig_mono in or remove it?
                                       'mods': None,
                                       'enzymes': None},
                   # TODO Make sure this isn't double counted
                   'assemble_ribosome_subunits': {'stoich': {'gtp_c': 1}
                                                  }}

ribosome_modifications = {'gtp_bound_30S_assembly_factor_phase1':
                          {'enzyme': 'Era_dim',
                           'stoich': {'gtp_c': 2,
                                      'h2o_c': 2,
                                      'h_c': -2,
                                      'pi_c': -2},
                           'num_mods': 1},

                          'RbfA_mono_assembly_factor_phase1':
                          {'enzyme': 'RbfA_mono',
                           'stoich': {},
                           'num_mods': 1},

                          'RimM_mono_assembly_factor_phase1':
                          {'enzyme': 'RimM_mono',
                           'stoich': {},
                           'num_mods': 1},

                          'Translation_initiation_factor_InfA':
                          {'enzyme': 'InfA_mono',
                           'stoich': {},
                           'num_mods': 1},

                          'Translation_initiation_factor_InfC':
                          {'enzyme': 'InfC_mono',
                           'stoich': {},
                           'num_mods': 1},

                          'Translation_gtp_initiation_factor_InfB':
                          {'enzyme': 'InfB_mono',
                           'stoich': {'gtp_c': 1,
                                      'h2o_c': 1,
                                      'h_c': -1,
                                      'pi_c': -1},
                           'num_mods': 1}}





def add_ribosome(me_model, verbose=True):
    ribosome_complex = ComplexData("ribosome", me_model)
    ribosome_components = ribosome_complex.stoichiometry

    for mod, components in rrna_modifications.items():
        rRNA_mod = ModificationData(mod, me_model)
        rRNA_mod.enzyme = components['machine']
        rRNA_mod.stoichiometry = components['metabolites']
        rRNA_mod.keff = 65.  # iOL uses 65. for all RNA mods
        if 'carriers' in components.keys():
            for carrier, stoich in components['carriers'].items():
                if stoich < 0:
                    rRNA_mod.enzyme += [carrier]
                rRNA_mod.stoichiometry[carrier] = stoich
        ribosome_complex.modifications[rRNA_mod.id] = 1

    mod_dict = ribosome_modifications
    for mod_id in mod_dict:
        mod_stoich = mod_dict[mod_id]['stoich']
        mod_enzyme = mod_dict[mod_id]['enzyme']
        num_mods = mod_dict[mod_id]['num_mods']
        mod = add_modification_data(me_model, mod_id, mod_stoich, mod_enzyme)
        ribosome_complex.modifications[mod.id] = -num_mods

    ribosome_assembly = ribosome_stoich
    for process in ribosome_assembly:
        for protein, amount in ribosome_assembly[process]['stoich'].items():
            ribosome_components[protein] += amount

    ribosome_complex.create_complex_formation(verbose=verbose)

# TODO: Double check the modifications here
rrna_modifications = {
                      # ---------16S Modifications---------------
                      'm2G_at_1207': {'machine': 'RsmC_mono',  # fixed
                                      'metabolites': {'amet_c': -1,
                                                      'ahcys_c': 1,
                                                      'h_c': 1}},
                      'm2G_at_1516': {'machine': None,
                                      # fixed, but still unknonw, NOT ybiN despite their ecocyc comments, the old 'MeT_16S_1516'
                                      'metabolites': {'amet_c': -1,
                                                      'ahcys_c': 1,
                                                      'h_c': 1}},
                      'm2G_at_966': {'machine': 'RsmD_mono',  # fixed
                                     'metabolites': {'amet_c': -1,
                                                     'ahcys_c': 1,
                                                     'h_c': 1}},
                      'm3U_at_1498': {'machine': 'YggJ_dim',  # fixed
                                      'metabolites': {'amet_c': -1,
                                                      'ahcys_c': 1,
                                                      'h_c': 1}},
                      'm4Cm_at_1402': {'machine': 'generic_16Sm4Cm1402',
                                       # fixed
                                       'metabolites': {'amet_c': -2,
                                                       'ahcys_c': 2,
                                                       'h_c': 2}},
                      'm5C_at_1407': {'machine': 'RsmF_mono',  # fixed
                                      'metabolites': {'amet_c': -1,
                                                      'ahcys_c': 1,
                                                      'h_c': 1}},
                      'm5C_at_967': {'machine': 'RsmB_mono',  # fixed
                                     'metabolites': {'amet_c': -1,
                                                     'ahcys_c': 1,
                                                     'h_c': 1}},
                      'm62A_at_1518': {'machine': 'KsgA_mono',  # fixed
                                       'metabolites': {'amet_c': -2,
                                                       'ahcys_c': 2,
                                                       'h_c': 2}},
                      'm62A_at_1519': {'machine': 'KsgA_mono',  # fixed
                                       'metabolites': {'amet_c': -2,
                                                       'ahcys_c': 2,
                                                       'h_c': 2}},
                      'm7G_at_527': {'machine': 'RsmG_mono',  # fixed
                                     'metabolites': {'amet_c': -1,
                                                     'ahcys_c': 1,
                                                     'h_c': 1}},
                      'Y_at_516': {'machine': 'RsuA_mono',  # fixed
                                   'metabolites': {}},

                      # ---------23S Modifications---------------
                      'Cm_at_2498': {'machine': 'RlmM_mono',
                                     'metabolites': {'amet_c': -1,
                                                     'ahcys_c': 1,
                                                     'h_c': 1}},
                      'D_at_2449': {'machine': None,
                                    # fixed, but still unknown, the old 'DU_23S_2449'
                                    'metabolites': {'h_c': -1,
                                                    'nadh_c': -1,
                                                    'nad_c': 1}},
                      'Gm_at_2251': {'machine': 'RlmB_dim',  # fixed
                                     'metabolites': {'amet_c': -1,
                                                     'ahcys_c': 1,
                                                     'h_c': 1}},
                      'm1G_at_745': {'machine': 'RrmA_dim_mod_2:zn2',
                                     # fixed
                                     'metabolites': {'amet_c': -1,
                                                     'ahcys_c': 1,
                                                     'h_c': 1}},
                      'm2A_at_2503': {'machine': 'RlmN_mono_mod_1:4fe4s',
                                      # fixed
                                      'metabolites': {'amet_c': -1,
                                                      'ahcys_c': 1,
                                                      'h_c': 1}},
                      'm2G_at_1835': {'machine': 'RlmG_mono',  # fixed
                                      'metabolites': {'amet_c': -1,
                                                      'ahcys_c': 1,
                                                      'h_c': 1}},
                      'm2G_at_2445': {'machine': 'RlmL_dim',  # fixed
                                      'metabolites': {'amet_c': -1,
                                                      'ahcys_c': 1,
                                                      'h_c': 1}},
                      'm5C_at_1962': {'machine': 'RlmI_dim',  # fixed
                                      'metabolites': {'amet_c': -1,
                                                      'ahcys_c': 1,
                                                      'h_c': 1}},
                      'm5U_at_1939': {'machine': 'RumA_mono_mod_1:4fe4s',
                                      # fixed
                                      'metabolites': {'amet_c': -1,
                                                      'ahcys_c': 1,
                                                      'h_c': 1}},
                      'm5U_at_747': {'machine': 'RumB_mono_mod_1:4fe4s',
                                     # fixed
                                     'metabolites': {'amet_c': -1,
                                                     'ahcys_c': 1,
                                                     'h_c': 1}},
                      'm6A_at_1618': {'machine': 'RlmF_mono',  # fixed
                                      'metabolites': {'amet_c': -1,
                                                      'ahcys_c': 1,
                                                      'h_c': 1}},
                      'm6A_at_2030': {'machine': None,
                                      # fixed, but still unknown, the old 'MeT_23S_2030'
                                      'metabolites': {'amet_c': -1,
                                                      'ahcys_c': 1,
                                                      'h_c': 1}},
                      'm7G_at_2069': {'machine': None,
                                      # fixed, but still unknonw, the old 'MeT_23S_2069'
                                      'metabolites': {'amet_c': -1,
                                                      'ahcys_c': 1,
                                                      'h_c': 1}},
                      'Um_at_2552': {'machine': 'RrmJ_mono',  # fixed
                                     'metabolites': {'amet_c': -1,
                                                     'ahcys_c': 1,
                                                     'h_c': 1}},
                      'Y_at_1911': {'machine': 'RluD_mono_mod_1:mg2',
                                    # fixed
                                    'metabolites': {}},
                      'Y_at_1915': {'machine': 'RluD_mono_mod_1:mg2',
                                    # fixed
                                    'metabolites': {}},
                      'm3Y_at_1915': {'machine': 'RlmH_dim',  # fixed
                                      'metabolites': {'amet_c': -1,
                                                      'ahcys_c': 1,
                                                      'h_c': 1}},
                      'Y_at_1917': {'machine': 'RluD_mono_mod_1:mg2',
                                    # fixed
                                    'metabolites': {}},
                      'Y_at_2457': {'machine': 'YmfC_mono',  # fixed
                                    'metabolites': {}},
                      'Y_at_2504': {'machine': 'RluC_mono',  # fixed
                                    'metabolites': {}},
                      'Y_at_2580': {'machine': 'RluC_mono',  # fixed
                                    'metabolites': {}},
                      'Y_at_2604': {'machine': 'YjbC_mono',  # fixed
                                    'metabolites': {}},
                      'Y_at_2605': {'machine': 'RluB_mono',  # fixed
                                    'metabolites': {}},
                      'Y_at_746': {'machine': 'RluA_mono',  # fixed
                                   'metabolites': {}},
                      'Y_at_955': {'machine': 'RluC_mono',  # fixed
                                   'metabolites': {}}}
