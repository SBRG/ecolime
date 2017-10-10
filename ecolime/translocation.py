from __future__ import print_function, absolute_import, division

from collections import defaultdict

from cobrame.core.processdata import PostTranslationData
from cobrame.core.reaction import PostTranslationReaction

pathway = {'sec': {'enzymes': {'SecB_tetra': {'length_dependent': True,
                                              'fixed_keff': False},
                               'SecA_MONOMER': {'length_dependent': True,
                                                'fixed_keff': False},
                               'Sec-CPLX': {'length_dependent': True,
                                            'fixed_keff': False}},
                   'keff': 4.,
                   'length_dependent_energy': True,
                   'stoichiometry': {'atp_c': -1./25., 'h2o_c': -1./25.,
                                     'adp_c': 1./25., 'pi_c': 1./25.,
                                     'h_c': 1./25.}},

           'tat': {'enzymes': {'TatBC_octa': {'length_dependent': False,
                                              'fixed_keff': False},
                               'TatA_MONOMER': {'length_dependent': False,
                                                'fixed_keff': False}},
                   'keff': 0.0125,
                   'length_dependent_energy': False,
                   'stoichiometry': {}},

           'tat_alt': {'enzymes': {'TatBC_octa': {'length_dependent': False,
                                                  'fixed_keff': False},
                                   'TatE_MONOMER': {'length_dependent': False,
                                                    'fixed_keff': False}},
                       'keff': 0.0125,
                       'length_dependent_energy': False,
                       'stoichiometry': {}},

           'bam': {'enzymes': {'SurA_MONOMER': {'length_dependent': False,
                                                'fixed_keff': False},
                               'Bam-CPLX': {'length_dependent': False,
                                            'fixed_keff': False}},
                   'keff': 0.027,
                   'length_dependent_energy': False,
                   'stoichiometry': {}},

           'lol': {'enzymes': {'LolCDE-CPLX': {'length_dependent': False,
                                               'fixed_keff': False},
                               'LolA_MONOMER': {'length_dependent': False,
                                                'fixed_keff': False},
                               'LolB_MONOMER': {'length_dependent': False,
                                                'fixed_keff': False}},
                   'keff': 0.9,
                   'length_dependent_energy': False,
                   'stoichiometry': {'atp_c': -1., 'h2o_c': -1.,
                                     'adp_c': 1., 'pi_c': 1., 'h_c': 1.}},

           'yidC': {'enzymes': {'SRP-CPLX': {'length_dependent': True,
                                             'fixed_keff': False},
                                'YidC_MONOMER': {'length_dependent': True,
                                                 'fixed_keff': False}},
                    'keff': 20.,
                    'length_dependent_energy': False,
                    'stoichiometry': {'gtp_c': -1., 'h2o_c': -1., 'gdp_c': 1.,
                                      'pi_c': 1., 'h_c': 1.}},

           'secA': {'enzymes': {'SecA_MONOMER': {'length_dependent': True,
                                                 'fixed_keff': False}},
                    'keff': 4.,
                    'length_dependent_energy': True,
                    'stoichiometry': {'atp_c': -1./3./25., 'h2o_c': -1./3./25.,
                                      'adp_c': 1./3./25., 'pi_c': 1./3./25.,
                                      'h_c': 1./3./25.}},

           'srp_yidC': {'enzymes': {'SRP-CPLX': {'length_dependent': True,
                                                 'fixed_keff': False},
                                    'YidC_MONOMER': {'length_dependent': True,
                                                     'fixed_keff': False},
                                    'Sec-CPLX': {'length_dependent': True,
                                                 'fixed_keff': False}},
                        'keff': 20.,
                        'length_dependent_energy': False,
                        'stoichiometry': {'gtp_c': -1., 'h2o_c': -1.,
                                          'gdp_c': 1., 'pi_c': 1., 'h_c': 1.}},

           'srp': {'enzymes': {'SRP-CPLX': {'length_dependent': True,
                                            'fixed_keff': False},
                               'FtsY_MONOMER': {'length_dependent': False,
                                                'fixed_keff': True},
                               'Sec-CPLX': {'length_dependent': True,
                                            'fixed_keff': False}},
                   'keff': 20.,
                   'length_dependent_energy': False,
                   'stoichiometry': {'gtp_c': -2., 'h2o_c': -2., 'gdp_c': 2.,
                                     'pi_c': 2., 'h_c': 2.}}
           }
abbreviation_to_pathway = {'s': 'sec_translocation',
                           't': ['tat_translocation', 'tat_translocation_alt'],
                           'b': 'bam_translocation',
                           'l': 'lol_translocation',
                           'y': 'yidC_translocation',
                           'a': 'secA_translocation',
                           'p': 'srp_yidC_translocation',
                           'r': 'srp_translocation'}

# Some proteins require different numbers of a complex in order to be
# translocated by a pathway
multipliers = {'YidC_MONOMER': {'b1855': 2., 'b3731': 2.},
               'TatE_MONOMER': {'b4072': 20.0, 'b1475': 5.166666666666667,
                                'b1474': 5.166666666666667, 'b2817': 22.0,
                                'b0973': 14.0, 'b0972': 14.0, 'b0997': 13.5,
                                'b3893': 15.5, 'b0152': 21.0, 'b3894': 15.5,
                                'b1589': 14.0, 'b1588': 15.0, 'b2938': 6.25,
                                'b1587': 14.0, 'b1872': 26.0, 'b2997': 14.0,
                                'b2996': 21.0, 'b2994': 14.0, 'b2435': 20.0,
                                'b0894': 26.0, 'b0123': 23.0, 'b2206': 27.0,
                                'b2205': 20.0},
               'TatA_MONOMER': {'b4072': 20.0, 'b1475': 5.166666666666667,
                                'b1474': 5.166666666666667, 'b2817': 22.0,
                                'b0973': 14.0, 'b0972': 14.0, 'b0997': 13.5,
                                'b3893': 15.5, 'b0152': 21.0, 'b3894': 15.5,
                                'b1589': 14.0, 'b1588': 15.0, 'b2938': 6.25,
                                'b1587': 14.0, 'b1872': 26.0, 'b2997': 14.0,
                                'b2996': 21.0, 'b2994': 14.0, 'b2435': 20.0,
                                'b0894': 26.0, 'b0123': 23.0, 'b2206': 27.0,
                                'b2205': 20.0}}

multipliers_protein_keys = defaultdict(dict)
for enzyme, value in multipliers.items():
    for bnum in value.keys():
        multipliers_protein_keys['protein_' + bnum][enzyme] = value[bnum]

mmol = 6.022e20  # number of molecules per mmol
nm2_per_m2 = 1e18  # used to convert nm^2 to m^2


def add_translocation_pathways(model, pathways_df, membrane_constraints=False):

    def add_translocation_data_and_reaction(model, pathways, preprocessed_id,
                                            processed_id, compartment,
                                            alt=False):

        suffix = '_alt' if alt else ''

        data = PostTranslationData('translocation_' + preprocessed_id + suffix,
                                   model, processed_id, preprocessed_id)

        data.translocation = pathways

        data.translocation_multipliers = \
            multipliers_protein_keys.get(preprocessed_id, {})

        # Add protein surface area constraint
        if membrane_constraints and compartment != 'Periplasm':
            mass = peptide_data.mass
            membrane_thickness = model.global_info['membrane_thickness']
            thickness = membrane_thickness[compartment]
            # Relationship uses protein molecular in kDa
            # Adds surface area constraint in units of m^2/mmol
            data.surface_area['SA_protein_' + compartment] = \
                (1.21 / thickness * 2.) * mass * mmol / nm2_per_m2

        rxn = PostTranslationReaction('translocation_' + peptide_data.id +
                                      suffix)
        rxn.posttranslation_data = data
        model.add_reaction(rxn)
        rxn.update()

    # loop through all translation data and add translocation rxns/surface area
    # constraints if they are membrane proteins
    for peptide_data in model.translation_data:

        # extract translocation info if peptide contained in complex
        # stoichiometry
        translocation_info = pathways_df[
            pathways_df.Protein.str.match(peptide_data.id)]
        # iterate if protein is not in a membrane complex
        if len(translocation_info) == 0:
            continue

        # Assign preprocessed and processed (translocated) peptide ids
        compartment = translocation_info.Protein_compartment.values[0]
        processed_id = 'protein_' + peptide_data.id + '_' + compartment
        preprocessed_id = 'protein_' + peptide_data.id

        # compile translocation pathways for each membrane protein
        pathways = set()
        pathways_alt = set()
        for abbrev in translocation_info.translocase_pathway.values[0]:
            pathway_name = abbreviation_to_pathway[abbrev]

            # The tat translocation pathway can use an alternate enzyme
            if type(pathway_name) == list:
                pathways.add(pathway_name[0])
                pathways_alt.add(pathway_name[1])
            else:
                pathways.add(pathway_name)
                pathways_alt.add(pathway_name)

        add_translocation_data_and_reaction(model, pathways, preprocessed_id,
                                            processed_id, compartment)

        # if there's an alternative pathway (tat) add this reaction as well
        if pathways != pathways_alt:
            add_translocation_data_and_reaction(model, pathways_alt,
                                                preprocessed_id, processed_id,
                                                compartment, alt=True)


lipoprotein_precursors = {'AcrA': 'b0463', 'AcrE': 'b3265', 'BamB': 'b2512',
                          'BamC': 'b2477', 'BamD': 'b2595', 'BamE': 'b2617',
                          'CusC': 'b0572', 'EG10544-MONOMER': 'b1677',
                          'EmtA': 'b1193', 'LolB': 'b1209', 'MetQ': 'b0197',
                          'MltA': 'b2813', 'MltB': 'b2701', 'MltC': 'b2963'}

lipid_modifications = ['pg120_p', 'pg141_p', 'pg140_p', 'pg181_p', 'pg161_p',
                       'pg160_p', 'pg180_p']


def add_lipoprotein_formation(model, compartment_dict,
                              membrane_constraints=False, update=True):

    # loop through all proteins which need lipid modifications (lipoproteins)
    for protein in lipoprotein_precursors.values():
        compartment = compartment_dict.get(protein)
        protein_met = model.metabolites.get_by_id('protein_' + protein)
        mass = protein_met.mass

        processed_id = 'protein_' + protein + '_lipoprotein_' + compartment
        preprocessed_id = 'protein_' + protein + '_' + compartment

        def add_lipoprotein_data_and_reaction(first_lipid, second_lipid):

            # Add PostTranslation Data, modifications and surface area
            data = PostTranslationData(reaction_prefix + '_' + second_lipid,
                                       model, processed_id, preprocessed_id)
            data.subreactions['mod_' + first_lipid] = 1
            data.subreactions['mod2_' + second_lipid + '_p'] = 1

            if membrane_constraints:
                thickness_dict = model.global_info['membrane_thickness']
                thickness = thickness_dict['Outer_Membrane']

                # From Liu et al. x2 for each to account for each leaflet
                protein_SA = 1.21 / thickness * 2 * mass * mmol / nm2_per_m2
                data.surface_area = {'SA_protein_' + compartment: -protein_SA,
                                     'SA_lipoprotein': 1. * mmol / nm2_per_m2}

            # Add Reaction to model and associated it with its data
            rxn = PostTranslationReaction(reaction_prefix + '_' + second_lipid)
            model.add_reaction(rxn)
            rxn.posttranslation_data = data
            if update:
                rxn.update()

        for mod in lipid_modifications:
            reaction_prefix = protein + '_lipid_modification_' + mod
            add_lipoprotein_data_and_reaction(mod, 'pg160')
            add_lipoprotein_data_and_reaction(mod, 'pe160')
