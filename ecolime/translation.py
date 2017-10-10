from __future__ import division, absolute_import, print_function

from cobrame import SubreactionData, Complex
from cobrame.util import dogma

# 1 machine + 1 atp + 1 aa + 1 h2o --> 1 machine-amp + 1 h + 1 ppi
# 1 machine-amp + 1 free tRNA --> 1 machine + 1 amp + 1 charged tRNA
special_trna_subreactions = {
    'sec_addition_at_UGA': {
        'enzymes': ['SelA_deca_mod_10:pydx5p',
                    'SelB_mono'],  # Selenocysteine loaders
        'stoich': {'h_c': 1, 'selnp_c': -1,
                   'pi_c': 1,
                   'generic_tRNA_UGA_cys__L_c': -1},
        'element_contribution': {'O': -1, 'Se': 1}}}

initiation_subreactions = {
    'Translation_initiation_factor_InfA':
        {'enzymes': 'InfA_mono',
         'stoich': {}},

    'Translation_initiation_factor_InfC':
        {'enzymes': 'InfC_mono',
         'stoich': {}},

    'Translation_gtp_initiation_factor_InfB':
        {'enzymes': 'InfB_mono',
         'stoich': {'gtp_c': -1,
                    'h2o_c': -1,
                    'h_c': 1,
                    'pi_c': 1,
                    'gdp_c': 1}},

    'fmet_addition_at_START':
        {'enzymes': ['InfB_mono',
                     'Fmt_mono_mod_mg2_mod_k'],
         # iOL had h_c:1 for fmet addition but this is not mass balanced
         'stoich': {'10fthf_c': -1, 'thf_c': 1,
                    # 'h_c': 1,
                    'generic_tRNA_START_met__L_c': -1},
         'element_contribution': {'C': 1, 'O': 1}}
   }

elongation_subreactions = {'FusA_mono_elongation': {'enzymes': ['FusA_mono'],
                                                    'stoich': {'gtp_c': -1,
                                                               'h2o_c': -1,
                                                               'h_c': 1,
                                                               'pi_c': 1,
                                                               'gdp_c': 1}},

                           'Tuf_gtp_regeneration': {'enzymes': ['Tsf_mono'],
                                                    'stoich': {}}}

# TODO Go through and double check elongation/termination ATP usage etc.
termination_subreactions = {'PrfA_mono_mediated_termination':
                            {'enzymes': ['PrfA_mono'],
                             'stoich': {}},

                            'PrfB_mono_mediated_termination':
                            {'enzymes': ['PrfB_mono'],
                             'stoich': {}},

                            'generic_RF_mediated_termination':
                            {'enzymes': ['generic_RF'],
                             'stoich': {}},

                            'N_terminal_methionine_cleavage':
                            {'enzymes': ['Map_mono_mod_2:fe2'],
                             'stoich': {'h2o_c': -1,
                                        'met__L_c': 1, 'h_c': 1},
                             'element_contribution': {'H': -10, 'O': -1,
                                                      'C': -5, 'N': -1,
                                                      'S': -1}},

                            'peptide_deformylase_processing':
                            {'enzymes': ['Def_mono_mod_1:fe2'],
                             'stoich': {'h2o_c': -1,
                                        'for_c': 1},
                             'element_contribution':
                                 {'H': 1, 'O': -1, 'C': -1}},

                            # This is a GTPS
                            'peptide_chain_release':
                            {'enzymes': ['PrfC_mono'],
                             'stoich': {'gtp_c': -1,
                                        'h2o_c': -1,
                                        'h_c': 1,
                                        'pi_c': 1,
                                        'gdp_c': 1}},

                            'ribosome_recycler':
                            {'enzymes': ['Rrf_mono'],
                             'stoich': {}},

                            'GroEL_dependent_folding':
                            {'enzymes': ['GroL_14', 'cisGroES_hepta',
                                         'transGroES_hepta'],
                             'stoich': {'atp_c': -7,
                                        'h2o_c': -7,
                                        'h_c': 7,
                                        'adp_c': 7,
                                        'pi_c': 7}},

                            # DnaK is correct
                            'DnaK_dependent_folding':
                            {'enzymes': ['DnaK_mono', 'DnaJ_dim_mod_4:zn2',
                                         'GrpE_dim'],
                             'stoich': {'atp_c': -1,
                                        'h2o_c': -1,
                                        'h_c': 1,
                                        'adp_c': 1,
                                        'pi_c': 1}}
                            }

# Subreaction for translation termination
translation_stop_dict = {'UAG': 'PrfA_mono',
                         'UGA': 'PrfB_mono',
                         'UAA': 'generic_RF'}

translation_start_codons = {"AUG", "GUG", "UUG", "AUU", "CUG"}

# Dictionary of frame shift mutations
frameshift_dict = {'b2891': '3033206:3034228,3034230:3034304'}

peptide_processing_subreactions = {"peptide_deformylase_processing",
                                   "peptide_chain_release",
                                   "ribosome_recycler"}


def add_translation_subreactions_to_model(me_model):
    # add general subreactions
    for rxn, info in elongation_subreactions.items():
        data = SubreactionData(rxn, me_model)
        data.enzyme = info['enzymes']
        data.stoichiometry = info['stoich']

    # add subreactions associated with termination and postprocessing
    for rxn, info in termination_subreactions.items():
        data = SubreactionData(rxn, me_model)
        data.enzyme = info['enzymes']
        data.stoichiometry = info['stoich']
        data._element_contribution = info.get('element_contribution', {})

    # add subreactions associated with translation initiation
    for rxn, info in initiation_subreactions.items():
        data = SubreactionData(rxn, me_model)
        data.enzyme = info['enzymes']
        data.stoichiometry = info['stoich']
        data._element_contribution = info.get('element_contribution', {})


def add_charged_trna_subreactions(me_model):
    # create subreaction for each codon. this will be used to model
    # the addition of charged tRNAs to the elongating peptide
    for codon in dogma.codon_table:
        if dogma.codon_table[codon] == '*':
            stop_codon = codon.replace('T', 'U')
            stop_enzyme = translation_stop_dict.get(stop_codon)
            me_model.add_metabolites([Complex(stop_enzyme)])

            subreaction_data = SubreactionData(
                stop_codon + '_' + stop_enzyme + '_mediated_termination',
                me_model)
            subreaction_data.enzyme = stop_enzyme
            subreaction_data.stoichiometry = {}
        else:
            full_aa = dogma.amino_acids[dogma.codon_table[codon]]
            amino_acid = full_aa.split('_')[0]
            subreaction_data = SubreactionData(
                amino_acid + '_addition_at_' + codon.replace('T', 'U'),
                me_model)
            trna = 'generic_tRNA_' + codon.replace('T', 'U') + '_' + full_aa
            subreaction_data.enzyme = 'generic_Tuf'  # Default AA loader enzyme

            # Accounts for GTP hydrolyzed by EF-TU and the ATP hydrolysis to
            # AMP required to add the amino acid to the tRNA
            subreaction_data.stoichiometry = {'gtp_c': -1, 'h2o_c': -2,
                                              'gdp_c': 1, 'h_c': 2, 'pi_c': 1,
                                              'ppi_c': 1, 'amp_c': 1,
                                              'atp_c': -1,
                                              trna: -1}

    # Add subreactions for start codon and selenocysteine
    for rxn, info in special_trna_subreactions.items():
        data = SubreactionData(rxn, me_model)
        data.enzyme = info['enzymes']
        data.stoichiometry = info['stoich']
        data._element_contribution = info.get('element_contribution', {})

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

# Codons are not unique to a tRNA
trna_to_codon = {'b2691': ['CGU', 'CGC', 'CGA'],
                 'b2692': ['CGU', 'CGC', 'CGA'],
                 'b2693': ['CGU', 'CGC', 'CGA'],
                 'b2694': ['CGU', 'CGC', 'CGA'], 'b2695': ['AGC', 'AGU'],
                 'b2590': ['GAA', 'GAG'], 'b0244': ['ACG'],
                 'b0883': ['UCC', 'UCU'], 'b1231': ['UAC', 'UAU'],
                 'b1230': ['UAC', 'UAU'], 'b4270': ['UUG'],
                 'b3853': ['GCA', 'GCG', 'GCU'], 'b3852': ['AUC', 'AUU'],
                 'b4164': ['GGC', 'GGU'], 'b4165': ['GGC', 'GGU'],
                 'b4370': ['CUG'], 'b4163': ['GGC', 'GGU'],
                 'b2189': ['CCC', 'CCU'], 'b1032': ['UCC', 'UCU'],
                 'b2348': ['AGG'], 'b3277': ['AUC', 'AUU'],
                 'b3276': ['GCA', 'GCG', 'GCU'], 'b3273': ['ACC', 'ACU'],
                 'b1986': ['AAC', 'AAU'], 'b0666': ['AUG'], 'b0665': ['CAA'],
                 'b0664': ['CAA'], 'b3969': ['GAA', 'GAG'], 'b4368': ['CUG'],
                 'b4369': ['CUG'], 'b2864': ['GGG'],
                 'b0971': ['UCU', 'UCA', 'UCG'], 'b0206': ['GAC', 'GAU'],
                 'b3757': ['GAA', 'GAG'], 'b2404': ['AAA', 'AAG'],
                 'b2401': ['GUA', 'GUG', 'GUU'],
                 'b2403': ['GUA', 'GUG', 'GUU'], 'b2396': ['GCC'],
                 'b3978': ['GGA'], 'b3979': ['ACC', 'ACU'], 'b2397': ['GCC'],
                 'b2402': ['GUA', 'GUG', 'GUU'],
                 'b3976': ['ACU', 'ACA', 'ACG'], 'b3977': ['UAC', 'UAU'],
                 'b3171': ['START'], 'b3174': ['CUC', 'CUU'], 'b1975': ['UCG'],
                 'b1977': ['AAC', 'AAU'], 'b3761': ['UGG'],
                 'b3760': ['GAC', 'GAU'], 'b4134': ['UUC', 'UUU'],
                 'b3545': ['CCG'], 'b0668': ['CAG'], 'b1984': ['AAC', 'AAU'],
                 'b0203': ['GCA', 'GCG', 'GCU'], 'b0202': ['AUC', 'AUU'],
                 'b1989': ['AAC', 'AAU'], 'b2967': ['UUC', 'UUU'],
                 'b3069': ['AUA'], 'b1909': ['UUG', 'UUA'], 'b0670': ['CAG'],
                 'b0672': ['CUG', 'CUA'], 'b0673': ['AUG'], 'b2816': ['START'],
                 'b2814': ['START'], 'b2815': ['START'],
                 'b0216': ['GAC', 'GAU'], 'b3797': ['CAC', 'CAU'],
                 'b4008': ['GAA', 'GAG'], 'b3798': ['CUG'],
                 'b3799': ['CCG', 'CCU', 'CCA'], 'b2652': ['AUA'],
                 'b1910': ['UGC', 'UGU'], 'b1911': ['GGC', 'GGU'],
                 'b0536': ['AGA'], 'b3658': ['UGA'], 'b1665': ['GUC', 'GUU'],
                 'b1666': ['GUC', 'GUU'], 'b3796': ['CGG'],
                 'b0744': ['GUA', 'GUG', 'GUU'], 'b0745': ['AAA', 'AAG'],
                 'b0746': ['GUA', 'GUG', 'GUU'], 'b0747': ['AAA', 'AAG'],
                 'b0743': ['AAA', 'AAG'], 'b0748': ['AAA', 'AAG'],
                 'b0749': ['AAA', 'AAG']}
