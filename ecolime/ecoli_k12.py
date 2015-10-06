__author__ = 'coltonlloyd'

divalent_list = ['ca2', 'mg2', 'mn2', 'cobalt2', 'ni2', 'cd2', 'zn2']
monovalent_list = ['k', 'na1']

generic_16s_rRNAs = ['b3851', 'b3968', 'b3756', 'b3278', 'b4007', 'b2591',
                     'b0201']
generic_23s_rRNAs = ['b3854', 'b3970', 'b3758', 'b3275', 'b4009', 'b2589',
                     'b0204']
generic_5s_rRNAs = ['b3855', 'b3971', 'b3759', 'b3274', 'b4010', 'b2588',
                    'b0205', 'b3272']

generic_RNase_list = ['RNase_T_dim_mod_4:mg2', 'RNase_BN_dim_mod_2:zn2',
                      'Rnd_mono_mod_5:mg2', 'Rnb_mono_mod_1:mg2',
                      'Rph_mono_mod_mg2']

Ribosome_modifications_phase1 = {'gtp_bound_30S_assembly_factor_phase1':
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
                                  'num_mods': 1}}

Ribosome_30s_proteins = {'RpsD_mono': 1,  # Phase 1 protein
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
                         'RpsA_mono': 1,  # Phase 2 proteins
                         'RpsB_mono': 1,
                         'RpsC_mono': 1,
                         'RpsJ_mono': 1,
                         'RpsN_mono': 1,
                         'RpsU_mono': 1,
                         'Sra_mono': 1}

Ribosome_50s_proteins = {'RplA_mono': 1,
                         'RplB_mono': 1,
                         'RplC_mono': 1,
                         'RplD_mono': 1,
                         'RplE_mono': 1,
                         'RplI_mono': 1,
                         'RplJ_mono': 1,  # TODO is the RplJ stoich right?
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
                         'rpL7/12_mod_1:acetyl': 2,
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
                         'RpmJ_mono': 1}

Ribosome_factors = {'Tig_mono': 1,
                    'InfA_mono': 1,
                    'InfC_mono': 1,
                    'InfB_mono': 1}  # adds to complex with gtp

# Other metabolites used in ribosome formation
Ribosome_metabolites = {'gtp_c': 1,
                        'mg2_c': 171}

Ribosome_rRNA_generics = {'generic_16s': 1,
                          'generic_23s': 1,
                          'generic_5s': 1}


def replace_divalent(c):
    for i in divalent_list:
        c = c.replace(i, "generic_divalent")
    return c

rRNA_containing = ['RNase_E_tetra_mod_2:zn2', 'RNase_P_cplx_mod_2:mg2',
                   'generic_RNase', 'RNase_m5', 'RNase_m16', 'RNase_m23',
                   'RNase_III_dim_mod_2:mg2', 'RNase_G_dim',
                   'RNase_T_dim_mod_4:mg2']

monocistronic = ['RNase_E_tetra_mod_2:zn2', 'RNase_P_cplx_mod_2:mg2',
                 'generic_RNase']

polycistronic_wout_rRNA = ['RNase_E_tetra_mod_2:zn2', 'RNase_P_cplx_mod_2:mg2',
                           'generic_RNase', 'RNase_III_dim', 'RNase_G_dim',
                           'RNase_T_dim_mod_4:mg2']

RNA_polymerase_components = {"b3295" : "rpoA",
                             "b3988" : "rpoC",
                             "b3987" : "rpoB"}