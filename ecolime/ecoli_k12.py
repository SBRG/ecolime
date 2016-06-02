
generic_RNase_list = ['RNase_T_dim_mod_4:mg2', 'RNase_BN_dim_mod_2:zn2',
                      'Rnd_mono_mod_5:mg2', 'Rnb_mono_mod_1:mg2',
                      'Rph_mono_mod_mg2']

generic_dict = {'generic_16Sm4Cm1402': ['RsmH_mono', 'RsmI_mono'],
                'generic_LYSINEaaRS': ['LysI_RS_dim',
                                       'LysII_RS_dim_mod_6:mg2'],
                'generic_Dus': ['DusA_mono', 'DusB_mono', 'DusC_mono'],
                'generic_RF': ['PrfA_mono', 'PrfB_mono'],
                'generic_Tuf': ['TufA_mono', 'TufB_mono'],
                'generic_RNase': ['RNase_T_dim_mod_4:mg2',
                                  'RNase_BN_dim_mod_2:zn2',
                                  'Rnd_mono_mod_5:mg2', 'Rnb_mono_mod_1:mg2',
                                  'Rph_mono_mod_mg2'],
                'generic_16s_rRNAs': ['RNA_b3851', 'RNA_b3968', 'RNA_b3756',
                                      'RNA_b3278', 'RNA_b4007', 'RNA_b2591',
                                      'RNA_b0201'],
                'generic_23s_rRNAs': ['RNA_b3854', 'RNA_b3970', 'RNA_b3758',
                                      'RNA_b3275', 'RNA_b4009', 'RNA_b2589',
                                      'RNA_b0204'],
                'generic_5s_rRNAs': ['RNA_b3855', 'RNA_b3971', 'RNA_b3759',
                                     'RNA_b3274', 'RNA_b4010', 'RNA_b2588',
                                     'RNA_b0205', 'RNA_b3272']}


excision_machinery = {
    'rRNA_containing': ['RNase_E_tetra_mod_2:zn2', 'RNase_P_cplx_mod_2:mg2',
                        'generic_RNase', 'RNase_m5', 'RNase_m16', 'RNase_m23',
                        'RNase_III_dim_mod_2:mg2', 'RNase_G_dim',
                        'RNase_T_dim_mod_4:mg2'],
    'monocistronic': ['RNase_E_tetra_mod_2:zn2', 'RNase_P_cplx_mod_2:mg2',
                      'generic_RNase'],
    'polycistronic_wout_rRNA': ['RNase_E_tetra_mod_2:zn2',
                                'RNase_P_cplx_mod_2:mg2', 'generic_RNase',
                                'RNase_III_dim', 'RNase_G_dim',
                                'RNase_T_dim_mod_4:mg2']}


polycistronic_wout_rRNA = ['RNase_E_tetra_mod_2:zn2', 'RNase_P_cplx_mod_2:mg2',
                           'generic_RNase', 'RNase_III_dim', 'RNase_G_dim',
                           'RNase_T_dim_mod_4:mg2']


RNA_polymerase_components = {"b3295": "rpoA",
                             "b3988": "rpoC",
                             "b3987": "rpoB"}

no_TU_list = ['b0024', 'b4586', 'b4690', 'b0533', 'b4588', 'b4589', 'b4590',
              'b0799', 'b4705', 'b4417', 'b0877', 'b0952', 'b1172', 'b4593',
              'b4594', 'b1181', 'b4672', 'b4699', 'b4674', 'b1413', 'b4601',
              'b4602', 'b4431', 'b4676', 'b4675', 'b4677', 'b1874', 'b4678',
              'b4436', 'b4667', 'b4437', 'b4668', 'b4679', 'b4604', 'b2218',
              'b4543', 'b4680', 'b2586', 'b4461', 'b4442', 'b4701', 'b4682',
              'b4443', 'b4683', 'b4684', 'b4446', 'b4665', 'b4664', 'b4666',
              'b3110', 'b3247', 'b3248', 'b4612', 'b4704', 'b3529', 'b4454',
              'b4557', 'b3705', 'b4456', 'b3865', 'b4686', 'b3933', 'b4620',
              'b4621', 'b4622', 'b4703', 'b4655']

RNA_degradosome = {'Eno_dim_mod_4:mg2': 1, 'Pnp_trim': 1,
                   'RNase_E_tetra_mod_2:zn2': 1, 'RhlB_dim': 1,
                   'Orn_dim_mod_2:mg2': 1}
