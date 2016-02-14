from minime import ComplexData, TranscribedGene, ModificationData
from minime.util.building import add_modification_data

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
                                                        'generic_16s_rRNAs': 1},
                                             'mods': {
                                                 'gtp_bound_30S_assembly_factor_phase1': 1},
                                             'enzymes': {'RbfA_mono',
                                                         'RimM_mono'}},
                   '30_S_assembly_2_(21S)': {'stoich': {'mg2_c': 60,
                                                        'RpsA_mono': 1,
                                                        'RpsB_mono': 1,
                                                        'RpsC_mono': 1,
                                                        'RpsJ_mono': 1, # TODO is the RplJ stoich right?
                                                        'RpsN_mono': 1,
                                                        'RpsU_mono': 1,
                                                        'Sra_mono': 1},
                                             'mods': None,
                                             'enzymes': None},
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
                                                  'rpL7/12_mod_1:acetyl': 2},
                                       'mods': None,
                                       'enzymes': None},
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
                   'assemble_ribosome_subunits': {'stoich': {'gtp_c': 1},
                                                  'mods': None, # maybe add InfB_gtp as modification not as separate process
                                                  'enzymes': {'InfB_mono',
                                                              'InfA_mono',
                                                              'InfC_mono'}}}

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
                           'stoich': {},
                           'num_mods': 1}}

# Subreaction for translation termination
translation_stop_dict = {'UAG': 'PrfA_mono',
                         'UGA': 'PrfB_mono',
                         'UAA': 'generic_RF'}

# 1 mRNA_nextRiboComplex + [1 b1211_assumedMonomer (pre-assigned) OR
#  1 b2891_assumedMonomer (pre-assigned) OR 1 generic_RF (pre-assigned)
#   depending on sequence ] + 1 b4375_assumedMonomer-gdp_FU.ID +
#  b0172_assumedMonomer (pre-assigned) + 1 gtp + 2 h2o -->
#  1 b3340_assumedMonomer-gdp +1 generic_Tuf-gdp + [1 b1211_assumedMonomer (not assigned) OR
#   1 b2891_assumedMonomer (not assigned) OR 1 generic_RF (not assigned)
#   depending on sequence] + 1  b4375_assumedMonomer-gdp + 3 h + 3 pi +
#  1 gdp+ 1 LAST tRNA uncharged + 1 peptide +1  b0172_assumedMonomer
#  (not assigned) + 1 ribosome (not assigned)
# TODO Go through and double check elongation/termination ATP usage etc.
translation_subreactions = {'PrfA_mono_mediated_termination':
                            {'enzyme': 'PrfA_mono',
                             'stoich': {}},

                            'PrfB_mono_mediated_termination':
                            {'enzyme': 'PrfB_mono',
                             'stoich': {}},

                            'generic_RF_mediated_termination':
                            {'enzyme': 'generic_RF',
                             'stoich': {}},

                            'N_terminal_methionine_cleavage':
                            {'enzyme': 'Map_mono_mod_2:fe2',
                             'stoich': {'h2o_c': -1,
                                        'met__L_c': 1}},

                            'peptide_deformylase_processing':
                            {'enzyme': 'Def_mono_mod_1:fe2',
                             'stoich': {'h2o_c': -1,
                                        'for_c': 1}},
                            # TODO look at this especially

                            'peptide_chain_release':
                            {'enzyme': 'PrfC_mono',
                             'stoich': {}},

                            'ribosome_recycler':
                            {'enzyme': 'Rrf_mono',
                             'stoich': {}},

                            # TODO below INCOMPLETE and this can be done better
                            # 7 adp and 7 mg2 are used to modify GroEL
                            'GroEL_dependent_folding':
                            {'enzyme': ['GroL_14', 'cisGroES_hepta',
                                        'transGroES_hepta'],
                             'stoich': {'atp_c': -7,
                                        'h2o_c': -7,
                                        'h_c': 7,
                                        'adp_c': 7,
                                        'pi_c': 7}},

                            'DnaK_dependent_folding':
                            {'enzyme': ['DnaK_mono', 'DnaJ_dim_mod_4:zn2',
                                        'GrpE_dim'],
                             'stoich': {'atp_c': -1,
                                        'h2o_c': -1,
                                        'h_c': 1,
                                        'adp_c': 1,
                                        'pi_c': 1}}
                            }

# Dictionary of frame shift mutations
frameshift_dict = {'b2891': '3033206:3034228,3034230:3034304'}

generic_rRNAs = {"generic_16s_rRNAs": ['b3851', 'b3968', 'b3756', 'b3278',
                                       'b4007', 'b2591', 'b0201'],
                 "generic_23s_rRNAs": ['b3854', 'b3970', 'b3758', 'b3275',
                                       'b4009', 'b2589', 'b0204'],
                 "generic_5s_rRNAs": ['b3855', 'b3971', 'b3759', 'b3274',
                                      'b4010', 'b2588', 'b0205', 'b3272']}


def add_ribosome(me_model, verbose=True):
    ribosome_complex = ComplexData("ribosome", me_model)
    ribosome_components = ribosome_complex.stoichiometry

    for mod, components in rrna_modifications.items():
        tRNA_mod = ModificationData(mod, me_model)
        tRNA_mod.enzyme = components['machine']
        tRNA_mod.stoichiometry = components['metabolites']
        if 'carriers' in components.keys():
            for carrier, stoich in components['carriers'].items():
                if stoich < 0:
                    tRNA_mod.enzyme += [carrier]
                tRNA_mod.stoichiometry[carrier] = stoich
        ribosome_complex.modifications[tRNA_mod.id] = 1

    for rRNA_type, generic_list in generic_rRNAs.items():
        for rRNA in generic_list:
            rRNA_id = 'RNA_' + rRNA
            me_model.add_metabolites([TranscribedGene(rRNA_type)])
            me_model.add_metabolites([TranscribedGene(rRNA_id)])
            new_rxn = cobra.Reaction("rRNA_" + rRNA + '_to_generic')
            me_model.add_reaction(new_rxn)
            new_rxn.reaction = rRNA_id + ' <=> ' + rRNA_type

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

# N terminal methionine cleaved
methionine_cleaved = ['b4154', 'b1109', 'b3908', 'b3417', 'b3940', 'b0344',
                      'b1263', 'b0026', 'b2697', 'b2114', 'b3559', 'b0680',
                      'b0033', 'b3341', 'b4120', 'b3742', 'b4254', 'b0932',
                      'b1923', 'b0918', 'b1062', 'b3635', 'b3771', 'b3774',
                      'b2133', 'b2566', 'b2831', 'b0115', 'b0273', 'b2020',
                      'b0051', 'b4201', 'b4013', 'b0092', 'b2499', 'b0063',
                      'b0015', 'b4024', 'b0074', 'b1224', 'b0121', 'b1960',
                      'b0903', 'b3213', 'b0072', 'b1094', 'b3390', 'b0388',
                      'b3172', 'b3731', 'b3926', 'b4143', 'b0438', 'b2750',
                      'b3287', 'b3225', 'b3340', 'b4147', 'b0170', 'b2779',
                      'b0954', 'b0757', 'b2153', 'b2904', 'b3430', 'b4172',
                      'b1712', 'b0014', 'b1718', 'b1207', 'b0783', 'b4219',
                      'b1778', 'b2518', 'b4162', 'b4245', 'b2564', 'b2926',
                      'b2743', 'b4041', 'b4226', 'b2501', 'b3932', 'b4177',
                      'b3642', 'b2780', 'b0171', 'b4244', 'b3167', 'b2699',
                      'b3700', 'b4375', 'b3985', 'b3983', 'b3986', 'b2606',
                      'b3984', 'b1716', 'b3185', 'b3637', 'b1089', 'b3636',
                      'b1717', 'b3297', 'b3342', 'b3298', 'b4202', 'b3316',
                      'b0023', 'b0169', 'b3314', 'b3296', 'b3303', 'b3306',
                      'b3230', 'b3649', 'b1863', 'b2942', 'b2620', 'b1324',
                      'b0008', 'b1261', 'b4393', 'b0779', 'b1004', 'b0422',
                      'b4362', 'b3780', 'b0930', 'b2890', 'b4129', 'b2962',
                      'b1937', 'b3938', 'b1187', 'b1584', 'b0623', 'b2097',
                      'b1779', 'b2927', 'b3870', 'b4062', 'b1661', 'b2479',
                      'b3247', 'b4207', 'b3775', 'b0439', 'b2297', 'b0902',
                      'b0116', 'b0888', 'b3162', 'b1241', 'b2525', 'b2913',
                      'b3201', 'b4391', 'b3556', 'b1823', 'b3781', 'b0058',
                      'b1333', 'b3725', 'b1095', 'b1092', 'b1589', 'b1003',
                      'b2925', 'b3770', 'b3732', 'b3739', 'b0185', 'b1912',
                      'b2414', 'b4042', 'b4384', 'b3639', 'b0812', 'b1203',
                      'b1073', 'b3198', 'b2303', 'b1612', 'b4153', 'b3610',
                      'b3229', 'b0369', 'b2994', 'b1237', 'b0889', 'b1658',
                      'b3753', 'b2597', 'b2665', 'b2898', 'b3764', 'b3165',
                      'b0605', 'b0168', 'b1882', 'b3495', 'b1288', 'b0028',
                      'b3844', 'b1236', 'b0507', 'b2787', 'b2231', 'b3699',
                      'b2717', 'b1175', 'b0782', 'b4243', 'b3982', 'b0727',
                      'b0114', 'b1638', 'b2502', 'b1304', 'b2312', 'b0523',
                      'b3302', 'b3305', 'b3307', 'b3311', 'b1656', 'b4059',
                      'b3089', 'b0729', 'b0452', 'b1637', 'b1482', 'b3339',
                      'b3980', 'b2801', 'b2329', 'b2580', 'b3831', 'b0776',
                      'b0778', 'b4019', 'b0088', 'b2908', 'b3822', 'b0237',
                      'b0160', 'b2708', 'b3729', 'b1734', 'b0090', 'b0312',
                      'b2763', 'b2762', 'b0931', 'b3359', 'b0895', 'b1468',
                      'b0060', 'b4179', 'b1854', 'b0179', 'b0674', 'b0085',
                      'b4293', 'b1702', 'b0907', 'b0467', 'b3288', 'b1538',
                      'b0529', 'b1924', 'b2563', 'b0240', 'b1276', 'b0337',
                      'b1415', 'b3829', 'b2156', 'b0142', 'b4039', 'b2905',
                      'b3962', 'b3924', 'b0109', 'b0073', 'b0071', 'b0785',
                      'b3176', 'b2521', 'b0036', 'b1967', 'b1896', 'b1461',
                      'b4386', 'b3856', 'b2903', 'b1849', 'b1082', 'b2286',
                      'b2283', 'b2594', 'b2186', 'b0147', 'b2976', 'b0828',
                      'b0368', 'b2764', 'b4371', 'b2675', 'b1412', 'b3081',
                      'b2687', 'b3437', 'b0714', 'b2785', 'b0059', 'b3317',
                      'b3309', 'b2514', 'b0684', 'b3308', 'b0755', 'b3065',
                      'b0884', 'b3839', 'b2417', 'b1817', 'b1591', 'b2235',
                      'b1076', 'b1118', 'b1396', 'b1519', 'b1919', 'b0334',
                      'b1182', 'b0420', 'b0331', 'b0459', 'b0935', 'b2747',
                      'b0484']

folding_dict = {
    'GroEL_dependent_folding': ['b0014', 'b0015', 'b0061', 'b0062', 'b0064',
                                'b0114', 'b0115', 'b0116', 'b0130', 'b0134',
                                'b0143', 'b0144', 'b0167', 'b0170', 'b0172',
                                'b0185', 'b0209', 'b0369', 'b0404', 'b0439',
                                'b0593', 'b0607', 'b0628', 'b0660', 'b0726',
                                'b0727', 'b0755', 'b0782', 'b0783', 'b0797',
                                'b0870', 'b0902', 'b0929', 'b0930', 'b0957',
                                'b1062', 'b1095', 'b1107', 'b1109', 'b1189',
                                'b1190', 'b1215', 'b1241', 'b1243', 'b1398',
                                'b1718', 'b1719', 'b1748', 'b1779', 'b1831',
                                'b2096', 'b2097', 'b2140', 'b2149', 'b2155',
                                'b2284', 'b2285', 'b2286', 'b2296', 'b2435',
                                'b2441', 'b2478', 'b2533', 'b2551', 'b2557',
                                'b2607', 'b2608', 'b2614', 'b2620', 'b2830',
                                'b2834', 'b2916', 'b2925', 'b2926', 'b2942',
                                'b3162', 'b3256', 'b3260', 'b3295', 'b3340',
                                'b3357', 'b3390', 'b3405', 'b3433', 'b3607',
                                'b3650', 'b3651', 'b3708', 'b3725', 'b3741',
                                'b3775', 'b3780', 'b3845', 'b3847', 'b3850',
                                'b3865', 'b3941', 'b3957', 'b3962', 'b3987',
                                'b3988', 'b3990', 'b3991', 'b3993', 'b3997',
                                'b4039', 'b4154', 'b4177', 'b4381', 'b4382'],
    'DnaK_dependent_folding': ['b2507', 'b2508', 'b2557', 'b2697', 'b2699',
                               'b2764', 'b2780', 'b2913', 'b2925', 'b2926',
                               'b2935', 'b3067', 'b3189', 'b3212', 'b3295',
                               'b3339', 'b3340', 'b3384', 'b3686', 'b3687',
                               'b3708', 'b3744', 'b3783', 'b3829', 'b3831',
                               'b3870', 'b3893', 'b3931', 'b3942', 'b3980',
                               'b3987', 'b3988', 'b4019', 'b4129', 'b4131',
                               'b4147', 'b4177', 'b4232', 'b4239', 'b4258',
                               'b4260', 'b4375', 'b4382', 'b4383', 'b4384',
                               'b4391', 'b0008', 'b0032', 'b0033', 'b0059',
                               'b0095', 'b0114', 'b0115', 'b0118', 'b0194',
                               'b0438', 'b0439', 'b0642', 'b0680', 'b0726',
                               'b0728', 'b0755', 'b0893', 'b0903', 'b0930',
                               'b0932', 'b1014', 'b1095', 'b1114', 'b1136',
                               'b1175', 'b1224', 'b1241', 'b1275', 'b1479',
                               'b1612', 'b1614', 'b1676', 'b1713', 'b1719',
                               'b1779', 'b1945', 'b2029', 'b2036', 'b2114',
                               'b2231', 'b2234', 'b2284', 'b2287', 'b2297',
                               'b2463']}
