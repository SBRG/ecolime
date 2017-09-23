from __future__ import print_function, absolute_import, division

from six import iteritems

import cobrame
from ecolime.corrections import correct_tRNA_modifications

amino_acid_tRNA_synthetase = {
  "cys__L_c": "CysS_mono_mod_1:zn2",
  "leu__L_c": "LeuS_mono",
  "lys__L_c": "generic_LYSINEaaRS",
  "asp__L_c": "Asp_RS_dim",
  "phe__L_c": "Phe_RS_tetra_mod_mg2",
  "his__L_c": "His_RS_dim_mod_4:mg2",
  "asn__L_c": "Asn_RS_dim",
  "pro__L_c": "Pro_RS_dim",
  "ala__L_c": "Ala_RS_tetra_mod_4:zn2",
  "ile__L_c": "IleS_mono_mod_2:zn2",
  "ser__L_c": "Ser_RS_dim_mod_mg2",
  "arg__L_c": "ArgS_mono_mod_mg2",
  "met__L_c": "Met_RS_dim_mod_2:zn2",
  "tyr__L_c": "Tyr_RS_dim",
  "glu__L_c": "GltX_mono_mod_mg2_mod_1:zn2",
  "thr__L_c": "Thr_RS_dim_mod_zn2",
  "val__L_c": "ValS_mono_mod_mg2",
  "gly_c": "Gly_RS_tetra",
  "trp__L_c": "Trp_RS_dim_mod_mg2",
  "gln__L_c": "GlnS_mono"
}

trna_modification = {'D_at_20A': {'machines': ['generic_Dus'],  
                                  'metabolites': {'nadph_c': -1,
                                                  'h_c': -1,
                                                  'nadp_c': 1}},

                     'D_at_20': {'machines': ['generic_Dus'],  
                                 'metabolites': {'nadph_c': -1,
                                                 'h_c': -1,
                                                 'nadp_c': 1}},

                     't6A_at_37': {'machines': ['YrdC_mono'],  
                                   'metabolites': {'hco3_c': -1,
                                                   'thr__L_c': -1,
                                                   'atp_c': -1,
                                                   'amp_c': 1,
                                                   'h_c': 1,
                                                   'h2o_c': 1,
                                                   'ppi_c': 1}},

                     'm7G_at_46': {'machines': ['YggH_mono'],  
                                   'metabolites': {'amet_c': -1,
                                                   'ahcys_c': 1,
                                                   'h_c': 1}},

                     'acp3U_at_47': {'machines': [],
                                     'metabolites': {'amet_c': -1,
                                                     '5mta_c': 1,
                                                     'h_c': 1}},

                     'm5U_at_54': {'machines': ['TrmA_mono'],  
                                   'metabolites': {'amet_c': -1,
                                                   'ahcys_c': 1,
                                                   'h_c': 1}},

                     'Y_at_55': {'machines': ['TruB_mono'],  
                                 'metabolites': {}},

                     'Y_at_65': {'machines': ['YqcB_mono'],  
                                 'metabolites': {}},

                     'D_at_17': {'machines': ['generic_Dus'],  
                                 'metabolites': {'nadph_c': -1,
                                                 'h_c': -1,
                                                 'nadp_c': 1}},

                     'cmo5U_at_34': {'machines': ['YecO_mono', 'YecP_mono'],
                                     'metabolites': {'amet_c': -1,
                                                     'chor_c': -2,
                                                     'ahcys_c': 1,
                                                     'h_c': 1,
                                                     'C10H8O5_c': 1,
                                                     'C9H9O4_c': 1}},

                     'D_at_16': {'machines': ['generic_Dus'],  
                                 'metabolites': {'nadph_c': -1,
                                                 'h_c': -1,
                                                 'nadp_c': 1}},

                     'Q_at_34': {
                     'machines': ['Tgt_hexa_mod_6:zn2', 'QueA_mono',
                                  'QueG_mono_mod_adocbl'],
                     'metabolites': {'preq1_c': -1,
                                     'amet_c': -1,
                                     'gua_c': 1,
                                     'ade_c': 1,
                                     'met__L_c': 1,
                                     'h_c': 2}},

                     'm2A_at_37': {'machines': [],
                                   'metabolites': {'amet_c': -1,
                                                   'ahcys_c': 1,
                                                   'h_c': 1}},

                     's4U_at_8': {'machines': ['ThiI_mono'],
                                  'carriers': {'trdrd_c': -1,
                                               'trdox_c': 1,
                                               'IscS_mod_2:pydx5p_mod_1:SH': -1,
                                               'IscS_mod_2:pydx5p': 1},
                                  'metabolites': {'atp_c': -1,
                                                  'amp_c': 1,
                                                  'ppi_c': 1,
                                                  'h_c': 1}},

                     'm6t6A_at_37': {'machines': ['YrdC_mono'],
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
                                                     'thf_c': 1}},

                     'Y_at_40': {'machines': ['TruA_dim'],  
                                 'metabolites': {}},

                     'Gm_at_18': {'machines': ['TrmH_dim'],  
                                  'metabolites': {'amet_c': -1,
                                                  'ahcys_c': 1,
                                                  'h_c': 1}},

                     'Um_at_32': {'machines': [],
                                  'metabolites': {'amet_c': -1,
                                                  'ahcys_c': 1,
                                                  'h_c': 1}},

                     'Y_at_38': {'machines': ['TruA_dim'],  
                                 'metabolites': {}},

                     'ac4C_at_34': {'machines': ['TmcA_mono'],  
                                    'metabolites': {'accoa_c': -1,
                                                    'coa_c': 1}},

                     'Y_at_39': {'machines': ['TruA_dim'],  
                                 'metabolites': {}},

                     # YhhP, YheLMN, YccK involved in sulfur transferase
                     # activity. TrmU catalyzes the addition of sulfur to
                     # uridine
                     'mnm5s2U_at_34': {'machines':
                                           ['TrmU_mono', 'YhhP_mono',
                                            'YheLMN_cplx', 'YccK_mono',
                                            'MnmEG_cplx_mod_fad_mod_2:k',
                                            'MnmC_mono_mod_fad'],
                                       'carriers': {
                                           'IscS_mod_2:pydx5p_mod_1:SH': -1,
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
                                   'metabolites': {'amet_c': -1,
                                                   'ahcys_c': 1,
                                                   'h_c': 1}},

                     'Cm_at_32': {'machines': ['TrmJ_dim'],
                                  'metabolites': {'amet_c': -1,
                                                  'ahcys_c': 1,
                                                  'h_c': 1}},

                     'ms2i6A_at_37': {'machines': ['MiaA_dim_mod_2:mg2',
                                                   'MiaB_mono_mod_1:4fe4s'],

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

                     'Y_at_32': {'machines': ['RluA_mono'],  
                                 'metabolites': {}},

                     'D_at_21': {'machines': ['generic_Dus'],  
                                 'metabolites': {'nadph_c': -1,
                                                 'h_c': -1,
                                                 'nadp_c': 1}},

                     'm1G_at_37': {'machines': ['TrmD_dim'],  
                                   'metabolites': {'amet_c': -1,
                                                   'ahcys_c': 1,
                                                   'h_c': 1}},

                     'Y_at_13': {'machines': ['TruD_mono'],  
                                 'metabolites': {}},

                     'k2C_at_34': {'machines': ['TilS_mono'],  
                                   'metabolites': {'atp_c': -1,
                                                   'lys__L_c': -1,
                                                   'ppi_c': 1,
                                                   'amp_c': 1,
                                                   'h_c': 2}},

                     'I_at_34': {'machines': ['TadA_dim_mod_2:zn2'],
                                 'metabolites': {'h2o_c': -1, 'h_c': -1,
                                                 'nh4_c': 1}},

                     'i6A_at_37': {'machines': ['MiaA_dim_mod_2:mg2'],  
                                   'metabolites': {'dmpp_c': -1,
                                                   'ppi_c': 1}},

                     'D_at_20_in_met_tRNA': {'machines': ['DusA_mono'],
                                             
                                             'metabolites': {'nadph_c': -1,
                                                             'h_c': -1,
                                                             'nadp_c': 1}},

                     'D_at_16_in_met_tRNA': {'machines': ['DusA_mono'],
                                             
                                             'metabolites': {'nadph_c': -1,
                                                             'h_c': -1,
                                                             'nadp_c': 1}},

                     'D_at_17_in_met_tRNA': {'machines': ['DusA_mono'],
                                             
                                             'metabolites': {'nadph_c': -1,
                                                             'h_c': -1,
                                                             'nadp_c': 1}},
                     'D_at_20A_in_met_tRNA': {'machines': ['DusA_mono'],
                                              
                                              'metabolites': {'nadph_c': -1,
                                                              'h_c': -1,
                                                              'nadp_c': 1}}
                     }

modification_info = {'D': {'elements': {'H': 2}, 'charge': 0},
                     'i6A': {'elements': {'C': 5, 'H': 8}, 'charge': 0},
                     'I': {'elements': {'N': -1, 'H': -1, 'O': 1},
                           'charge': 0},
                     'k2C': {'elements': {'O': 1, 'N': 2, 'H': 12, 'C': 6},
                             'charge': 0},
                     'Y': {'elements': {}, 'charge': 0},
                     'm1G': {'elements': {'H': 2, 'C': 1}, 'charge': 0},
                     'ms2i6A': {'elements': {'C': 6, 'H': 10, 'S': 1},
                                'charge': 0},
                     'Cm': {'elements': {'H': 2, 'C': 1}, 'charge': 0},
                     'Um': {'elements': {'H': 2, 'C': 1}, 'charge': 0},
                     'm6A': {'elements': {'H': 2, 'C': 1}, 'charge': 0},
                     'mnm5s2U': {'elements': {'C': 2, 'H': 5, 'N': 1, 'O': -1,
                                              'S': 1},
                                 'charge': 0},
                     'ac4C': {'elements': {'H': 2, 'C': 2, 'O': 1},
                              'charge': 0},
                     'Gm': {'elements': {'H': 2, 'C': 1}, 'charge': 0},
                     'mnm5U': {'elements': {'C': 2, 'H': 5, 'N': 1, 'O': -1,
                                            'S': 1},
                               'charge': 0},
                     's2C': {'elements': {'O': -1, 'S': 1}, 'charge': 0},
                     'm6t6A': {'elements': {'C': 6, 'O': 4, 'N': 1, 'H': 9},
                               'charge': 0},
                     's4U': {'elements': {'O': -1, 'S': 1}, 'charge': 0},
                     'm2A': {'elements': {'H': 2, 'C': 1}, 'charge': 0},
                     'Q': {'elements': {'C': 7, 'O': 2, 'H': 11}, 'charge': 1},
                     'cmo5U': {'elements': {'C': 3, 'O': 2, 'H': 4},
                               'charge': 0},
                     'm5U': {'elements': {'C': 1, 'H': 2}, 'charge': 0},
                     'acp3U': {'elements': {'C': 4, 'H': 7, 'N': 1, 'O': 2},
                               'charge': 0},
                     'm7G': {'elements': {'C': 1, 'H': 2}, 'charge': 0},
                     't6A': {'elements': {'C': 5, 'N': 1, 'O': 4, 'H': 6},
                             'charge': 0}
                     }


def add_tRNA_modification_procedures(model):

    modifications = trna_modification.copy()
    modifications = correct_tRNA_modifications(modifications)

    for mod, components in iteritems(modifications):
        tRNA_mod = cobrame.ModificationData(mod, model)
        tRNA_mod.enzyme = components['machines']
        tRNA_mod.stoichiometry = components['metabolites']
        tRNA_mod.keff = 65.  # iOL uses 65 for all tRNA mods
        if 'carriers' in components.keys():
            for carrier, stoich in components['carriers'].items():
                if stoich < 0:
                    tRNA_mod.enzyme += [carrier]
                tRNA_mod.stoichiometry[carrier] = stoich

        # Add element contribution from modification to tRNA
        tRNA_mod._element_contribution = \
            modification_info[mod.split('_')[0]]['elements']

    return modifications
