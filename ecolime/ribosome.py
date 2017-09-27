from __future__ import print_function, absolute_import, division

from six import iteritems

from cobrame import ComplexData, SubreactionData
from ecolime.corrections import correct_rrna_modifications


# Dictionary of {formation_step:[{metabolite:stoichiometry}]}
# Positive for reactants negative for products (complex formation convention)

# Leaving out Tig_mono trigger factor for now. It's a chaperone not part of
# formation

###### From Teddy's notes in original ME code #######
# What I'm trying to accomplish here: (old, but still useful)
# 1) 1 30Sp + 1 generic_16S -> 1 rib_30
# 2) 1 50Sp + 1 generic_23S + 1 generic_5S -> 1 rib_50
# 3) 1 rib_30 + 1 rib_50 -> 1 rib_70
# 4) 1 rib_70 + 1 b0884_assumedMonomer + 1 b1718_uniprotComplex  -->
#       1 rib_30_IF1_IF3 + 1 rib_50
# 5) 1 b3168_assumedMonomer_gtp + 1 rib_30_IF1_IF3 --> 1 rib_30_ini


# Mod:
#   Era_dim (assembly factor) + 2 gtp +2 h20->
#   30S_assembly_factor_gtp_hydrolying_assembly_phase_1_gtp ->
#           2 gtp + 2 h + 2 pi


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
                                                        'generic_16s_rRNAs':
                                                            1}},
                   '30_S_assembly_2_(21S)': {'stoich': {'RpsA_mono': 1,
                                                        'RpsB_mono': 1,
                                                        'RpsC_mono': 1,
                                                        'RpsJ_mono': 1,
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
                   '50_S_assembly_2': {'stoich': {'RplF_mono': 1,
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
                                                  'Tig_mono': 1},
                                       'mods': None,
                                       'enzymes': None},
                   # TODO Make sure this isn't double counted
                   'assemble_ribosome_subunits': {'stoich': {'gtp_c': 1}
                                                  }}

ribosome_subreactions = {'gtp_bound_30S_assembly_factor_phase1':
                         {'enzyme': 'Era_dim',
                          'stoich': {'gtp_c': -2,
                                     'h2o_c': -2,
                                     'h_c': 2,
                                     'pi_c': 2,
                                     'gdp_c': 2},
                          'num_mods': 1},

                         'RbfA_mono_assembly_factor_phase1':
                         {'enzyme': 'RbfA_mono',
                          'stoich': {},
                          'num_mods': 1},

                         'RimM_mono_assembly_factor_phase1':
                         {'enzyme': 'RimM_mono',
                          'stoich': {},
                          'num_mods': 1}
                         }


def add_ribosome(me_model, verbose=True):
    ribosome_complex = ComplexData("ribosome", me_model)
    ribosome_components = ribosome_complex.stoichiometry

    rrna_mods = correct_rrna_modifications(rrna_modifications)
    for mod, components in iteritems(rrna_mods):
        rrna_mod = SubreactionData(mod, me_model)
        rrna_mod.enzyme = components['machine']
        rrna_mod.stoichiometry = components['metabolites']
        rrna_mod.keff = 65.  # iOL uses 65. for all RNA mods

        # Add element contribution from modification to rRNA
        rrna_mod._element_contribution = \
            modification_info[mod.split('_')[0]]['elements']

        if 'carriers' in components.keys():
            for carrier, stoich in iteritems(components['carriers']):
                if stoich < 0:
                    rrna_mod.enzyme += [carrier]
                rrna_mod.stoichiometry[carrier] = stoich
        ribosome_complex.subreactions[rrna_mod.id] = 1

    subreaction_dict = ribosome_subreactions
    for subreaction_id in subreaction_dict:
        # get subreaction info
        subreaction_stoich = subreaction_dict[subreaction_id]['stoich']
        subreaction_enzyme = subreaction_dict[subreaction_id]['enzyme']
        num_subreactions = subreaction_dict[subreaction_id]['num_mods']

        # add subreaction to model
        subreaction = SubreactionData(subreaction_id, me_model)
        subreaction.stoichiometry = subreaction_stoich
        subreaction.enzyme = subreaction_enzyme

        # account for subreactions in complex data
        ribosome_complex.subreactions[subreaction.id] = num_subreactions

    # Ribosomes in iOL1650 contain 171 mg2 ions
    ribosome_complex.subreactions['mod_mg2_c'] = 171.
    ribosome_assembly = ribosome_stoich
    for process in ribosome_assembly:
        for protein, amount in iteritems(ribosome_assembly[process]['stoich']):
            ribosome_components[protein] += amount

    ribosome_complex.create_complex_formation(verbose=verbose)

rrna_modifications = {
                      # ---------16S Modifications---------------
                      'm2G_at_1207': {'machine': 'RsmC_mono',  
                                      'metabolites': {'amet_c': -1,
                                                      'ahcys_c': 1,
                                                      'h_c': 1}},
                      'm2G_at_1516': {'machine': None,
                                      'metabolites': {'amet_c': -1,
                                                      'ahcys_c': 1,
                                                      'h_c': 1}},
                      'm2G_at_966': {'machine': 'RsmD_mono',  
                                     'metabolites': {'amet_c': -1,
                                                     'ahcys_c': 1,
                                                     'h_c': 1}},
                      'm3U_at_1498': {'machine': 'YggJ_dim',  
                                      'metabolites': {'amet_c': -1,
                                                      'ahcys_c': 1,
                                                      'h_c': 1}},
                      'm4Cm_at_1402': {'machine': 'generic_16Sm4Cm1402',
                                       'metabolites': {'amet_c': -2,
                                                       'ahcys_c': 2,
                                                       'h_c': 2}},
                      'm5C_at_1407': {'machine': 'RsmF_mono',  
                                      'metabolites': {'amet_c': -1,
                                                      'ahcys_c': 1,
                                                      'h_c': 1}},
                      'm5C_at_967': {'machine': 'RsmB_mono',  
                                     'metabolites': {'amet_c': -1,
                                                     'ahcys_c': 1,
                                                     'h_c': 1}},
                      'm62A_at_1518': {'machine': 'KsgA_mono',  
                                       'metabolites': {'amet_c': -2,
                                                       'ahcys_c': 2,
                                                       'h_c': 2}},
                      'm62A_at_1519': {'machine': 'KsgA_mono',  
                                       'metabolites': {'amet_c': -2,
                                                       'ahcys_c': 2,
                                                       'h_c': 2}},
                      'm7G_at_527': {'machine': 'RsmG_mono',  
                                     'metabolites': {'amet_c': -1,
                                                     'ahcys_c': 1,
                                                     'h_c': 1}},
                      'Y_at_516': {'machine': 'RsuA_mono',  
                                   'metabolites': {}},

                      # ---------23S Modifications---------------
                      'Cm_at_2498': {'machine': 'RlmM_mono',
                                     'metabolites': {'amet_c': -1,
                                                     'ahcys_c': 1,
                                                     'h_c': 1}},
                      'D_at_2449': {'machine': None,
                                    'metabolites': {'h_c': -1,
                                                    'nadh_c': -1,
                                                    'nad_c': 1}},
                      'Gm_at_2251': {'machine': 'RlmB_dim',  
                                     'metabolites': {'amet_c': -1,
                                                     'ahcys_c': 1,
                                                     'h_c': 1}},
                      'm1G_at_745': {'machine': 'RrmA_dim_mod_2:zn2',
                                     'metabolites': {'amet_c': -1,
                                                     'ahcys_c': 1,
                                                     'h_c': 1}},
                      'm2A_at_2503': {'machine': 'RlmN_mono_mod_1:4fe4s',
                                      'metabolites': {'amet_c': -1,
                                                      'ahcys_c': 1,
                                                      'h_c': 1}},
                      'm2G_at_1835': {'machine': 'RlmG_mono',  
                                      'metabolites': {'amet_c': -1,
                                                      'ahcys_c': 1,
                                                      'h_c': 1}},
                      'm2G_at_2445': {'machine': 'RlmL_dim',  
                                      'metabolites': {'amet_c': -1,
                                                      'ahcys_c': 1,
                                                      'h_c': 1}},
                      'm5C_at_1962': {'machine': 'RlmI_dim',  
                                      'metabolites': {'amet_c': -1,
                                                      'ahcys_c': 1,
                                                      'h_c': 1}},
                      'm5U_at_1939': {'machine': 'RumA_mono_mod_1:4fe4s',
                                      'metabolites': {'amet_c': -1,
                                                      'ahcys_c': 1,
                                                      'h_c': 1}},
                      'm5U_at_747': {'machine': 'RumB_mono_mod_1:4fe4s',
                                     'metabolites': {'amet_c': -1,
                                                     'ahcys_c': 1,
                                                     'h_c': 1}},
                      'm6A_at_1618': {'machine': 'RlmF_mono',  
                                      'metabolites': {'amet_c': -1,
                                                      'ahcys_c': 1,
                                                      'h_c': 1}},
                      'm6A_at_2030': {'machine': None,
                                      'metabolites': {'amet_c': -1,
                                                      'ahcys_c': 1,
                                                      'h_c': 1}},
                      'm7G_at_2069': {'machine': None,
                                      'metabolites': {'amet_c': -1,
                                                      'ahcys_c': 1,
                                                      'h_c': 1}},
                      'Um_at_2552': {'machine': 'RrmJ_mono',  
                                     'metabolites': {'amet_c': -1,
                                                     'ahcys_c': 1,
                                                     'h_c': 1}},
                      'Y_at_1911': {'machine': 'RluD_mono_mod_1:mg2',
                                    'metabolites': {}},
                      'Y_at_1915': {'machine': 'RluD_mono_mod_1:mg2',
                                    'metabolites': {}},
                      'm3Y_at_1915': {'machine': 'RlmH_dim',  
                                      'metabolites': {'amet_c': -1,
                                                      'ahcys_c': 1,
                                                      'h_c': 1}},
                      'Y_at_1917': {'machine': 'RluD_mono_mod_1:mg2',
                                    'metabolites': {}},
                      'Y_at_2457': {'machine': 'YmfC_mono',  
                                    'metabolites': {}},
                      'Y_at_2504': {'machine': 'RluC_mono',  
                                    'metabolites': {}},
                      'Y_at_2580': {'machine': 'RluC_mono',  
                                    'metabolites': {}},
                      'Y_at_2604': {'machine': 'YjbC_mono',  
                                    'metabolites': {}},
                      'Y_at_2605': {'machine': 'RluB_mono',  
                                    'metabolites': {}},
                      'Y_at_746': {'machine': 'RluA_mono',  
                                   'metabolites': {}},
                      'Y_at_955': {'machine': 'RluC_mono',  
                                   'metabolites': {}}}

modification_info = {'Y': {'elements': {}, 'charge': 0},
                     'm3Y': {'elements': {'C': 1, 'H': 2}, 'charge': 0},
                     'Um': {'elements': {'H': 2, 'C': 1}, 'charge': 0},
                     'm7G': {'elements': {'C': 1, 'H': 2}, 'charge': 0},
                     'm6A': {'elements': {'C': 1, 'H': 2}, 'charge': 0},
                     'm5U': {'elements': {'C': 1, 'H': 2}, 'charge': 0},
                     'm5C': {'elements': {'C': 1, 'H': 2}, 'charge': 0},
                     'm2G': {'elements': {'C': 1, 'H': 2}, 'charge': 0},
                     'm2A': {'elements': {'C': 1, 'H': 2}, 'charge': 0},
                     'm1G': {'elements': {'C': 1, 'H': 2}, 'charge': 0},
                     'Gm': {'elements': {'C': 1, 'H': 2}, 'charge': 0},
                     'D': {'elements': {'H': 2}, 'charge': 0},
                     'Cm': {'elements': {'C': 1, 'H': 2}, 'charge': 0},
                     'm62A': {'elements': {'C': 2, 'H': 4}, 'charge': 0},
                     'm4Cm': {'elements': {'C': 2, 'H': 4}, 'charge': 0},
                     'm3U': {'elements': {'C': 1, 'H': 2}, 'charge': 0}
                     }
