
trna_modification = {'D_at_20A': {'machines': ['generic_Dus'],  # fixed
                                  'metabolites': {'nadph_c': -1,
                                                  'h_c': -1,
                                                  'nadp_c': 1}},

                     'D_at_20': {'machines': ['generic_Dus'],  # fixed
                                 'metabolites': {'nadph_c': -1,
                                                 'h_c': -1,
                                                 'nadp_c': 1}},

                     't6A_at_37': {'machines': ['YrdC_mono'],  # fixed
                                   'metabolites': {'hco3_c': -1,
                                                   'thr__L_c': -1,
                                                   'atp_c': -1,
                                                   'amp_c': 1,
                                                   'h_c': 1,
                                                   'h2o_c': 1,
                                                   'ppi_c': 1}},

                     'm7G_at_46': {'machines': ['YggH_mono'],  # fixed
                                   'metabolites': {'amet_c': -1,
                                                   'ahcys_c': 1,
                                                   'h_c': 1}},

                     'acp3U_at_47': {'machines': [],
                                     # fixed, still unknown, previously called 'AcpT_tRNA_pos_47_acp3U'
                                     'metabolites': {'amet_c': -1,
                                                     '5mta_c': 1,
                                                     'h_c': 1}},

                     'm5U_at_54': {'machines': ['TrmA_mono'],  # fixed
                                   'metabolites': {'amet_c': -1,
                                                   'ahcys_c': 1,
                                                   'h_c': 1}},

                     'Y_at_55': {'machines': ['TruB_mono'],  # fixed
                                 'metabolites': {}},

                     'Y_at_65': {'machines': ['YqcB_mono'],  # fixed
                                 'metabolites': {}},

                     'D_at_17': {'machines': ['generic_Dus'],  # fixed
                                 'metabolites': {'nadph_c': -1,
                                                 'h_c': -1,
                                                 'nadp_c': 1}},

                     'cmo5U_at_34': {'machines': ['YecO_mono', 'YecP_mono'],
                                     # fixed, also includes an unknown factor 'HyL_tRNA_pos_34_ho5U'
                                     'metabolites': {'amet_c': -1,
                                                     'chor_c': -2,
                                                     'ahcys_c': 1,
                                                     'h_c': 1,
                                                     'C10H8O5_c': 1,
                                                     # new metabolite warning!
                                                     'C9H9O4_c': 1}},
                     # new metabolite warning!

                     'D_at_16': {'machines': ['generic_Dus'],  # fixed
                                 'metabolites': {'nadph_c': -1,
                                                 'h_c': -1,
                                                 'nadp_c': 1}},

                     'Q_at_34': {
                     'machines': ['Tgt_hexa_mod_6:zn2', 'QueA_mono',
                                  'QueG_mono_mod_adocbl'],
                     # fixed and gapped the unknown factor
                     'metabolites': {'preq1_c': -1,  # yes
                                     'amet_c': -1,  # yes
                                     'gua_c': 1,  # yes
                                     'ade_c': 1,  # yes
                                     'met__L_c': 1,  # yes
                                     'h_c': 2}},  # yes

                     'm2A_at_37': {'machines': [],
                                   # fixed, but still unknown, previously called 'MeT_tRNA_pos_37_m2A'
                                   'metabolites': {'amet_c': -1,
                                                   'ahcys_c': 1,
                                                   'h_c': 1}},

                     's4U_at_8': {'machines': ['ThiI_mono'],
                                  # still need names for carriers
                                  'carriers': {'trdrd_c': -1,
                                               'trdox_c': 1,
                                               'IscS_mod_2:pydx5p_mod_1:SH': -1,
                                               'IscS_mod_2:pydx5p': 1},
                                  'metabolites': {'atp_c': -1,
                                                  'amp_c': 1,
                                                  'ppi_c': 1,
                                                  'h_c': 1}},

                     'm6t6A_at_37': {'machines': ['YrdC_mono'],
                                     # fixed, but this reaction also has a methylation by an unknown factor 'MeT_tRNA_pos_37_m6t6A', although the evidence for this should be checked more extensively
                                     'metabolites': {'amet_c': -1,
                                                     'atp_c': -1,
                                                     'hco3_c': -1,
                                                     'thr__L_c': -1,
                                                     'ahcys_c': 1,
                                                     'amp_c': 1,
                                                     'h_c': 2,
                                                     'h2o_c': 1,
                                                     'ppi_c': 1}},

                     's2C_at_32': {'machines': ['YdaO_mono'],
                                   # still need names for carriers
                                   'carriers': {'trdrd_c': -1,
                                                'trdox_c': 1,
                                                'IscS_mod_2:pydx5p_mod_1:SH': -1,
                                                'IscS_mod_2:pydx5p': 1},
                                   'metabolites': {'atp_c': -1,
                                                   'amp_c': 1,
                                                   'ppi_c': 1,
                                                   'h_c': 1}},

                     'mnm5U_at_34': {'machines': ['MnmEG_cplx_mod_fad_mod_2:k',
                                                  'MnmC_mono_mod_fad'],
                                     # fixed, but check the role of fad w/ mnmG and TrmCand fadh2 (not sure why fadh2 was in the original reaction).. maybe we can ask harish about this one. # CHECK that MnmEG is actually formed once Teddy adds the complex to the database.
                                     'metabolites': {'gtp_c': -1,
                                                     'h2o_c': -1,
                                                     '5fthf_c': -1,
                                                     'gly_c': -1,
                                                     'amet_c': -1,
                                                     'ahcys_c': 1,
                                                     'h_c': 3,
                                                     'gdp_c': 1,
                                                     'glx_c': 1,
                                                     'pi_c': 1,
                                                     # was missing this!
                                                     'thf_c': 1}},

                     'Y_at_40': {'machines': ['TruA_dim'],  # fixed
                                 'metabolites': {}},

                     'Gm_at_18': {'machines': ['TrmH_dim'],  # fixed
                                  'metabolites': {'amet_c': -1,
                                                  'ahcys_c': 1,
                                                  'h_c': 1}},

                     'Um_at_32': {'machines': [],
                                  # fixed, but still unknown, previously called 'MeT_tRNA_pos_32_Um'
                                  'metabolites': {'amet_c': -1,
                                                  'ahcys_c': 1,
                                                  'h_c': 1}},

                     'Y_at_38': {'machines': ['TruA_dim'],  # fixed
                                 'metabolites': {}},

                     'ac4C_at_34': {'machines': ['TmcA_mono'],  # fixed
                                    'metabolites': {'accoa_c': -1,
                                                    'coa_c': 1}},

                     'Y_at_39': {'machines': ['TruA_dim'],  # fixed
                                 'metabolites': {}},

                     'mnm5s2U_at_34': {'machines': ['TrmU_mono', 'YhhP_mono',
                                                    'YheLMN_cplx', 'YccK_mono',
                                                    'MnmEG_cplx_mod_fad_mod_2:k',
                                  'MnmC_mono_mod_fad'],
                     # check that MnmEG_cplex and YheLMN_cplx are formed. also check role of fadh2 (not sure why fadh2 was in the original reaction).
                     'carriers': {'IscS_mod_2:pydx5p_mod_1:SH': -1,
                                  'trdrd_c': -1,
                                  'IscS_mod_2:pydx5p': 1,
                                  'trdox_c': 1},
                     'metabolites': {'atp_c': -1,
                                     'gtp_c': -1,
                                     'h2o_c': -1,
                                     '5fthf_c': -1,
                                     'gly_c': -1,
                                     'amet_c': -1,
                                     'gdp_c': 1,
                                     'pi_c': 1,
                                     'h_c': 4,
                                     'thf_c': 1,
                                     'glx_c': 1,
                                     'ahcys_c': 1,
                                     'amp_c': 1,
                                     'ppi_c': 1}},

                     'm6A_at_37': {'machines': ['YfiC_mono'],
                                   # fixed, new gene
                                   'metabolites': {'amet_c': -1,
                                                   'ahcys_c': 1,
                                                   'h_c': 1}},

                     'Cm_at_32': {'machines': ['TrmJ_dim'],
                                  # fixed, previously this was called 'MeT_tRNA_pos_32_Cm'
                                  'metabolites': {'amet_c': -1,
                                                  'ahcys_c': 1,
                                                  'h_c': 1}},

                     'ms2i6A_at_37': {'machines': ['MiaA_dim_mod_2:mg2',
                                                   'MiaB_mono_mod_1:4fe4s'],
                                      # still need names for carriers
                                      'carriers': {
                                      'IscS_mod_2:pydx5p_mod_1:SH': -1,
                                      'IscS_mod_2:pydx5p': 1,
                                      'fldrd_c': -1,
                                      'fldox_c': 1, },
                                      'metabolites': {'dmpp_c': -1,
                                                      'amet_c': -2,
                                                      'ppi_c': 1,
                                                      'ahcys_c': 1,
                                                      'h_c': 2,
                                                      'met__L_c': 1,
                                                      'dad__5_c': 1,
                                                      }},

                     'Y_at_32': {'machines': ['RluA_mono'],  # fixed
                                 'metabolites': {}},

                     'D_at_21': {'machines': ['generic_Dus'],  # fixed
                                 'metabolites': {'nadph_c': -1,
                                                 'h_c': -1,
                                                 'nadp_c': 1}},

                     'm1G_at_37': {'machines': ['TrmD_dim'],  # fixed
                                   'metabolites': {'amet_c': -1,
                                                   'ahcys_c': 1,
                                                   'h_c': 1}},

                     'Y_at_13': {'machines': ['TruD_mono'],  # fixed
                                 'metabolites': {}},

                     'k2C_at_34': {'machines': ['TilS_mono'],  # fixed
                                   'metabolites': {'atp_c': -1,
                                                   'lys__L_c': -1,
                                                   'ppi_c': 1,
                                                   'amp_c': 1,
                                                   'h_c': 2}},

                     'I_at_34': {'machines': ['TadA_dim_mod_2:zn2'],
                                 # fixed (added an h to both sides, so nh4 is released instead of nh3)
                                 'metabolites': {'h2o_c': -1, 'h_c': -1,
                                                 'nh4_c': 1}},

                     'i6A_at_37': {'machines': ['MiaA_dim_mod_2:mg2'],  # fixed
                                   'metabolites': {'dmpp_c': -1,
                                                   'ppi_c': 1}},

                     'D_at_20_in_met_tRNA': {'machines': ['DusA_mono'],
                                             # fixed
                                             'metabolites': {'nadph_c': -1,
                                                             'h_c': -1,
                                                             'nadp_c': 1}},

                     'D_at_16_in_met_tRNA': {'machines': ['DusA_mono'],
                                             # fixed
                                             'metabolites': {'nadph_c': -1,
                                                             'h_c': -1,
                                                             'nadp_c': 1}},

                     'D_at_17_in_met_tRNA': {'machines': ['DusA_mono'],
                                             # fixed
                                             'metabolites': {'nadph_c': -1,
                                                             'h_c': -1,
                                                             'nadp_c': 1}},
                     'D_at_20A_in_met_tRNA': {'machines': ['DusA_mono'],
                                              # fixed
                                              'metabolites': {'nadph_c': -1,
                                                              'h_c': -1,
                                                              'nadp_c': 1}}
                     }


def get_tRNA_modification_procedures():

    mod = trna_modification.copy()

    # flavodoxin fix based off of doi:10.1016/j.febslet.2005.05.047
    # "Two types of flavodoxins exist in E. coli. Flavodoxin
    # I is encoded by fldA and is constitutively expressed in
    # E. coli whereas flavodoxin II is encoded by fldB and is induced
    # by oxidative stress [18]. Flavodoxin I has been shown
    # to be an essential gene in E. coli. Flavodoxin II cannot replace
    # flavodoxin I"

    correct_mod = mod['ms2i6A_at_37']['carriers']
    correct_mod["FLAVODOXIN1-MONOMER"] = correct_mod.pop('fldrd_c')
    correct_mod["FLAVODOXIN1-MONOMER_mod_Oxidized"] = correct_mod.pop('fldox_c')

    # iron sulfur clusters do sulfur transferase
    # no indication that that a separate redox metabolite is needed
    correct_mod = mod['s2C_at_32']['carriers']
    correct_mod.pop('trdrd_c')
    correct_mod.pop('trdox_c')

    # also not needed for this one. PMID 10753862
    correct_mod = mod['s4U_at_8']['carriers']
    correct_mod.pop('trdrd_c')
    correct_mod.pop('trdox_c')

    # no reference for this, but it's free anyway
    correct_mod = mod['mnm5s2U_at_34']['carriers']
    correct_mod.pop('trdrd_c')
    correct_mod.pop('trdox_c')

    return mod
