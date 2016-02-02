from minime import ComplexData, TranscribedGene
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
                                                        'RpsJ_mono': 1,
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

translation_subreactions = {'PrfA_mono_mediated_termination':
                            {'enzyme': 'PrfA_mono',
                             'stoich': {},
                             'num_subreactions': 1},
                            'PrfB_mono_mediated_termination':
                            {'enzyme': 'PrfB_mono',
                             'stoich': {},
                             'num_subreactions': 1},
                            'generic_RF_mediated_termination':
                            {'enzyme': 'generic_RF',
                             'stoich': {},
                             'num_subreactions': 1}}

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
    ribosome_modifications = ribosome_complex.modifications

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
        ribosome_modifications[mod.id] = -num_mods

    ribosome_assembly = ribosome_stoich
    for process in ribosome_assembly:
        for protein, amount in ribosome_assembly[process]['stoich'].items():
            ribosome_components[protein] += amount

    ribosome_complex.create_complex_formation(verbose=verbose)
