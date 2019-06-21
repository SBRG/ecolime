from __future__ import print_function, division, absolute_import

import math
import re

from . import flat_files
import cobrame
from cobrame.core.processdata import SubreactionData, PostTranslationData
from cobrame.core.component import ProcessedProtein, TranslatedGene
from cobrame.core.reaction import PostTranslationReaction


folding_subreactions = {
    'folding_KJE_1': {'enzymes': ['DnaJ_dim_mod_4:zn2',
                                  'DnaK_mono_bound_to_atp'],
                      'stoichiometry': {'h2o_c': -1,
                                        'h_c': 1,
                                        'pi_c': 1},
                      'keff': 0.5},
    'folding_KJE_4': {'enzymes': ['DnaJ_dim_mod_4:zn2',
                                  'DnaK_mono_bound_to_atp'],
                      'stoichiometry': {'h2o_c': -1,
                                        'h_c': 1,
                                        'pi_c': 1},
                      'keff': 0.5},

    'folding_KJE_2': {'enzymes': ['DnaK_mono_bound_to_atp'],
                      'stoichiometry': {'atp_c': -1,
                                        'adp_c': 1},
                      'keff': 0.04},

    'folding_KJE_GrpE': {'enzymes': ['GrpE_dim'],
                      'stoichiometry': {},
                      'keff': 0.2},

    'folding_KJE_3': {'enzymes': ['DnaK_mono_bound_to_atp'],
                      'stoichiometry': {'atp_c': -1,
                                        'adp_c': 1},
                      'keff': 0.04},

    'folding_KJE_5': {'enzymes': ['DnaK_mono_bound_to_atp'],
                      'stoichiometry': {'atp_c': -1,
                                        'adp_c': 1},
                      'keff': 0.04},

    'folding_KJE_6': {'enzymes': ['DnaK_mono_bound_to_atp'],
                      'stoichiometry': {'atp_c': -1,
                                        'adp_c': 1},
                      'keff': 0.04},

    'folding_GroEL_ES_1': {'enzymes':['[GroL]14[GroS]7_cis_with_7_adp_and_7_mg2','transGroES_hepta'],
                           'stoichiometry': {},
                           'keff': 50.0},
    'folding_GroEL_ES_4': {'enzymes':['[GroL]14[GroS]7_cis_with_7_adp_and_7_mg2','transGroES_hepta'],
                           'stoichiometry': {},
                           'keff': 50.0},

    'folding_GroEL_ES_2': {'enzymes':['[GroL]14[GroS]7_cis_with_7_adp_and_7_mg2'],
                           'stoichiometry': {'atp_c': -7,
                                             'h2o_c': -7,
                                             'adp_c': 7,
                                             'h_c': 7,
                                             'pi_c': 7},
                           'keff': 0.12},

    'folding_GroEL_ES_3': {'enzymes':['[GroL]14[GroS]7_cis_with_7_adp_and_7_mg2'],
                           'stoichiometry': {'atp_c': -7,
                                             'h2o_c': -7,
                                             'adp_c': 7,
                                             'h_c': 7,
                                             'pi_c': 7},
                           'keff': 0.12},

    'folding_GroEL_ES_5': {'enzymes':['[GroL]14[GroS]7_cis_with_7_adp_and_7_mg2'],
                           'stoichiometry': {'atp_c': -7,
                                             'h2o_c': -7,
                                             'adp_c': 7,
                                             'h_c': 7,
                                             'pi_c': 7},
                           'keff': 0.12},

    'folding_GroEL_ES_6': {'enzymes':['[GroL]14[GroS]7_cis_with_7_adp_and_7_mg2'],
                           'stoichiometry': {'atp_c': -7,
                                             'h2o_c': -7,
                                             'adp_c': 7,
                                             'h_c': 7,
                                             'pi_c': 7},
                           'keff': 0.12},

    #'folding_spontaneous': {'enzymes': None, 'stoichiometry': {}}
    'folding_spontaneous': {'enzymes': ['Lon'],
                            'stoichiometry': {'atp_c': -2,
                                              'h2o_c': -2,
                                              'adp_c': 2,
                                              'h_c': 2,
                                              'pi_c': 2},
                            'keff': 1.225},
    'folding_Lon': {'enzymes': ['Lon'],
                    'stoichiometry': {'atp_c': -2,
                                      'h2o_c': -2,
                                      'adp_c': 2,
                                      'h_c': 2,
                                      'pi_c': 2},
                    'keff': 1.225}

    }

genes_to_use_dill = ["b2411", "b4035", "b1709", "b2926", "b0115", "b1676",
                     "b2750", "b0033", "b3041", "b2414", "b1723", "b0915",
                     "b0179", "b2416", "b3236", "b2508", "b3187", "b0529",
                     "b2153", "b0133", "b2478", "b0008", "b0171", "b0945",
                     "b0521", "b3803", "b0825", "b2373", "b3199", "b0697",
                     "b0696", "b1252", "b2224", "b0590", "b2241", "b3290",
                     "b1677", "b0463", "b0572", "b2963", "b3469", "b3734",
                     "b1817", "b3731", "b3732", "b3736", "b3733", "b3735",
                     "b0484", "b3729", "b0149", "b0086", "b3730", "b0084",
                     "b3189", "b3176", "b4169", "b2134", "b2435", "b2027",
                     "b4161", "b3946", "b2458", "b2150", "b3566", "b1395",
                     "b3124", "b3125", "b1397", "b3567", "b1900", "b3600",
                     "b2169", "b2914", "b4087", "b2417", "b0509", "b0331",
                     "b0207", "b2415", "b4268", "b2542", "b3749", "b3751",
                     "b2297", "b1198", "b3748", "b1393", "b3752", "b2149",
                     "b2221", "b0350", "b0351", "b2341", "b2342", "b3541",
                     "b3540", "b2306", "b3455", "b3117", "b1014", "b2957",
                     "b2440", "b0365", "b4260", "b2779", "b1854", "b4395",
                     "b3956", "b0727", "b1611", "b3061", "b3479", "b3480",
                     "b2311", "b3450", "b2579", "b3779", "b0936", "b4094",
                     "b0383", "b3608", "b2529", "b0933", "b2751", "b1090",
                     "b2197", "b2201", "b3255", "b0617", "b1094", "b3198",
                     "b2708", "b3052", "b2711", "b0523", "b0323", "b4244",
                     "b1207", "b0910", "b2747", "b0029", "b4039", "b2515",
                     "b0421", "b1858", "b1210", "b3368", "b3805", "b1993",
                     "b1991", "b0009", "b0765", "b0783", "b0103", "b3639",
                     "b2564", "b0109", "b1740", "b0066", "b3992", "b2103",
                     "b0414", "b0415", "b0185", "b1093", "b3918", "b0809",
                     "b2329", "b1693", "b2601", "b3390", "b0388", "b4024",
                     "b3959", "b0864", "b2818", "b0243", "b0242", "b2677",
                     "b0166", "b0655", "b2421", "b2599", "b1704", "b2024",
                     "b2021", "b2020", "b2019", "b3771", "b0078", "b0074",
                     "b3925", "b3222", "b3223", "b4467", "b2418", "b0469",
                     "b3650", "b2066", "b4381", "b1098", "b2065", "b0125",
                     "b2498", "b2428", "b0615", "b0616", "b1325", "b2303",
                     "b2874", "b3035", "b4287", "b3201", "b0199", "b0588",
                     "b4485", "b1483", "b2129", "b0829", "b3200", "b0592",
                     "b0158", "b4227", "b0197", "b0574", "b0308", "b3288",
                     "b1713", "b3591", "b0893", "b1646", "b3265", "b2954",
                     "b0650", "b2530", "b0614", "b2744", "b0474", "b3993",
                     "b3503", "b2817", "b2261", "b0980", "b1922", "b3305",
                     "b4200", "b4203", "b3341", "b3985", "b3306", "b3230",
                     "b3983", "b1717", "b3321", "b3986", "b3297", "b3299",
                     "b3231", "b3342", "b3185", "b3310", "b3298", "b3301",
                     "b3307", "b3313", "b3165", "b3294", "b2609", "b3304",
                     "b3311", "b2606", "b4202", "b1716", "b3315", "b0023",
                     "b3318", "b3065", "b3309", "b3637", "b3312", "b3302",
                     "b1480", "b1089", "b3636", "b3703", "b3984", "b3317",
                     "b0911", "b3320", "b3314", "b3319", "b3296", "b3303",
                     "b2567", "b1084", "b3704", "b3643", "b3164", "b3202",
                     "b2741", "b4171", "b3166", "b4143", "b3167", "b3247",
                     "b1211", "b4142", "b0014", "b2891", "b3783", "b2573",
                     "b0172", "b2526", "b1086", "b0170", "b3181", "b2207",
                     "b2614", "b3706", "b3169", "b3295", "b0416", "b3987",
                     "b3988", "b3067", "b3461", "b0436", "b3649", "b3590",
                     "b0884", "b2594", "b3168", "b1718", "b2946", "b3651",
                     "b3343", "b3146", "b1135", "b4180", "b3465", "b2532",
                     "b4022", "b2785", "b0082", "b1269", "b2517", "b2140",
                     "b3741", "b3287", "b2566", "b0914", "b0194", "b3389",
                     "b0002", "b2913", "b1779", "b4388", "b3919", "b1264",
                     "b2610", "b0098", "b3464"]

groel_targets = ['b3340', 'b3988', 'b1241', 'b1719', 'b2286', 'b2557', 'b3650',
                 'b0726', 'b0115', 'b0114', 'b0014', 'b3987', 'b1243', 'b3741',
                 'b4154', 'b2155']


def add_chaperone_subreactions(model):
    for subreaction, values in folding_subreactions.items():
        data = SubreactionData(subreaction, model)
        data.enzyme = values['enzymes']
        data.stoichiometry = values['stoichiometry']
        #if subreaction != 'folding_spontaneous':
        data.keff = values['keff']  # in 1/s


def add_chaperone_network(model):
    pd = 0.25
    pg = 0.50

    # First remove folding subreactions from translation reactions
    for data in model.translation_data:
        for subreaction in list(data.subreactions.keys()):
            if 'folding' in subreaction:
                data.subreactions.pop(subreaction)

    dill_df = flat_files.get_dill_keq_df()
    oobatake_df = flat_files.get_oobatake_keq_df()
    k_folding_df = flat_files.get_folding_rates_df()
    aggregation_propensity_df = flat_files.get_aggregation_popensity_df()

    # Get set of all proteins contained in the model
    model_protein_set = set()
    for protein in model.metabolites.query(re.compile('^protein_b[0-9]')):
        if not isinstance(protein, TranslatedGene):
            continue
        model_protein_set.add(protein.id.split('_')[1])

    # For each protein, create
    for protein_bnum in list(model_protein_set):
        protein = model.metabolites.get_by_id('protein_'+protein_bnum)
        folded_met = ProcessedProtein(protein.id + '_folded', protein.id)
        model.add_metabolites([folded_met])

        for suffix in ['_intermediate', '_KJE_folding_intermediate1',
                       '_KJE_folding_intermediate2',
                       '_GroEL_ES_folding_intermediate1',
                       '_GroEL_ES_folding_intermediate2']:
            intermediate_met = ProcessedProtein(protein.id + suffix,
                                                protein.id)
            model.add_metabolites([intermediate_met])

        # TODO clean this up
        try:
            if protein_bnum in genes_to_use_dill:
                keq_folding = dill_df.T[protein_bnum].to_dict()
            else:
                keq_folding = oobatake_df.T[protein_bnum].to_dict()

            k_folding = k_folding_df.T[protein_bnum].to_dict()

            propensity = aggregation_propensity_df['propensity'][protein_bnum]
            if propensity == 0.:
                propensity = .01

        except:
            propensity = aggregation_propensity_df['propensity'].median()
            keq_folding = oobatake_df.median().to_dict()
            k_folding = k_folding_df.median().to_dict()

        for folding in folding_subreactions:
            if '_GrpE' in folding:
                continue
            folding_id = 'folding_' + protein.id + '_' + folding

            data = PostTranslationData(folding_id, model,
                                       protein.id + '_folded', protein.id)
            data.folding = True

            data.subreactions[folding] = 1
            if 'KJE' in folding and folding.split('_')[-1] not in ['1', '4']:
                data.subreactions['folding_KJE_GrpE'] = 1
            data.folding_mechanism = folding
            data.aggregation_propensity = propensity
            data.k_folding = k_folding
            data.keq_folding = keq_folding
            data.biomass_type = 'prosthetic_group_biomass'
            if folding.split('_')[-1] in []:
                data.intermediate_folding_id = ''
            else:
                data.intermediate_folding_id = None

            if protein.id in groel_targets or \
                    protein.formula_weight / 1000 < 60:
                data.size_or_target_scaling_factor = 2/7

            if 'KJE' in folding:
                kind = 'KJE'
                multiplier = pd
            elif 'GroEL_ES' in folding:
                kind = 'GroEL_ES'
                multiplier = pg

            if folding == 'folding_Lon':
                data.folding_reactant_metabolite = protein.id + '_intermediate'

                data.coupling_expression_id = 'folding_protease'
                data.dilution_multiplier = 1 - pd

            elif folding == 'folding_spontaneous':
                data.folding_reactant_metabolite = protein.id
                data.folding_product_metabolite = folded_met.id

                data.coupling_expression_id = \
                    'folding_spontaneous_w_degradation'

            if folding.endswith('_1'):
                data.folding_reactant_metabolite = protein.id
                data.folding_product_metabolite = \
                    '%s_%s_folding_intermediate1' % (protein.id, kind)
                data.coupling_expression_id = 'noncoupling'

            elif folding.endswith('_2'):
                data.folding_reactant_constraint = \
                    '%s_folding_%s_constraint1' % (protein.id, kind)
                data.folding_product_metabolite = folded_met.id
                data.folding_reactant_metabolite = \
                    '%s_%s_folding_intermediate1' % (protein.id, kind)

                data.dilution_multiplier = multiplier
                data.coupling_expression_id = 'folding_chaperone'

            elif folding.endswith('_3'):
                data.folding_product_metabolite = protein.id + '_intermediate'
                data.folding_product_constraint = \
                    '%s_folding_%s_constraint1' % (protein.id, kind)
                data.folding_reactant_metabolite = \
                    '%s_%s_folding_intermediate1' % (protein.id, kind)

                data.coupling_expression_id = 'noncoupling'

            elif folding.endswith('_4'):
                data.folding_reactant_metabolite = protein.id + '_intermediate'
                data.folding_product_metabolite =\
                    '%s_%s_folding_intermediate2' % (protein.id, kind)

                data.coupling_expression_id = 'noncoupling'

            elif folding.endswith('_5'):
                data.folding_product_metabolite = folded_met.id
                data.folding_reactant_metabolite = \
                    '%s_%s_folding_intermediate2' % (protein.id, kind)
                data.folding_reactant_constraint = \
                    '%s_folding_%s_constraint2' % (protein.id, kind)

                data.dilution_multiplier = 1 - multiplier
                data.coupling_expression_id = 'folding_chaperone'

            elif folding.endswith('_6'):
                data.folding_product_metabolite = protein.id + '_intermediate'
                data.folding_product_constraint = \
                    '%s_folding_%s_constraint2' % (protein.id, kind)
                data.folding_reactant_metabolite = \
                    '%s_%s_folding_intermediate2' % (protein.id, kind)

                data.coupling_expression_id = 'noncoupling'

            rxn = PostTranslationReaction(folding_id)
            model.add_reaction(rxn)
            rxn.posttranslation_data = data

            rxn.update()


def change_temperature(model, temperature):
    """
    Updates temperature dependent model parameters

    :param model: ME-model
    :param temperature: int (in degC)
    :return:
    """
    int_temp = int(temperature)
    str_temp = str(temperature)
    model.global_info['temperature'] = int_temp

    def get_temperature_dependent_keff(keff_0, T):
        T_kelvin = T + 273.15
        R = 1.9872  # in cal / mol / k
        T0 = 273.15 + 37.  # in k
        Ea = R * T0 * (30.5 - math.log(keff_0))

        return math.e ** (30.5 - Ea / R / T_kelvin)

    for rxn in model.reactions:
        if hasattr(rxn, 'keff'):
            new_keff = get_temperature_dependent_keff(rxn.keff, temperature)
            rxn.keff = new_keff
    for data in model.process_data:
        if hasattr(data, 'keff'):
            new_keff = get_temperature_dependent_keff(data.keff, temperature)
            data.keff = new_keff
        if hasattr(data, 'synthetase_keff'):
            new_keff = get_temperature_dependent_keff(data.synthetase_keff,
                                                      temperature)
            data.synthetase_keff = new_keff

    model.update()


def add_foldme_module(me):
    for r in me.subreaction_data:
        if not hasattr(r, 'keff'):
            r.keff = 65.

    rxn = cobrame.MEReaction('cisGroES_to_transGroES')
    me.add_reaction(rxn)
    rxn.add_metabolites({'cisGroES_hepta': -1, 'transGroES_hepta': 1})
    rxn.lower_bound = -1000.

    mod = cobrame.SubreactionData('mod_adp_c', me)
    mod.stoichiometry = {'adp_c': -1.0}

    data = cobrame.ComplexData('[GroL]14[GroS]7_cis_with_7_adp_and_7_mg2', me)
    data.stoichiometry = {'GroL_14': 1., 'cisGroES_hepta': 1.}
    data.subreactions = {'mod_adp_c': 7., 'mod_mg2_c': 7.}
    data.create_complex_formation()

    # modification to DnaK_mono
    mod = cobrame.SubreactionData('mod_atp_c', me)
    mod.stoichiometry = {'atp_c': -1.0}
    data = cobrame.ComplexData('DnaK_mono_bound_to_atp', me)
    data.stoichiometry = {'DnaK_mono': 1.}
    data.subreactions = {'mod_atp_c': 1.}
    data.create_complex_formation()

    # create Lon complex
    data = cobrame.ComplexData('Lon', me)
    data.stoichiometry = {'protein_b0439': 6}  # TODO protein actually b0439
    data.subreactions = {'mod_mg2_c': 6.}
    data.create_complex_formation()

    # add folding/cycling coupling constraint
    for met in me.metabolites.query(re.compile('RNA_b[0-9]')):
        protein_bnum = met.id.replace('RNA_', '')
        constraint_id_list = []
        constraint_id_list.append(
            'protein_' + protein_bnum + '_folding_KJE_constraint1')
        constraint_id_list.append(
            'protein_' + protein_bnum + '_folding_KJE_constraint2')
        constraint_id_list.append(
            'protein_' + protein_bnum + '_folding_GroEL_ES_constraint1')
        constraint_id_list.append(
            'protein_' + protein_bnum + '_folding_GroEL_ES_constraint2')
        for constraint_id in constraint_id_list:
            try:
                folding_constraint = me.metabolites.get_by_id(constraint_id)
            except KeyError:
                folding_constraint = cobrame.Constraint(constraint_id)
                me.add_metabolites([folding_constraint])

    # constructing chaperone network
    add_chaperone_subreactions(me)
    add_chaperone_network(me)

    # Update stoichiometries of complexes to the folded version of the proteins,
    # if they are not membrane/periplasm proteins
    for data in list(me.complex_data):
        for key in data.stoichiometry:
            if key == 'protein_dummy':
                continue
            if key.endswith('Membrane') or key.endswith('Periplasm'):
                continue
            if key.startswith('protein_') and not key.endswith('_folded'):
                value = data.stoichiometry.pop(key)
                data.stoichiometry[key + '_folded'] = value

    for data in me.complex_data:
        data.formation.update()

