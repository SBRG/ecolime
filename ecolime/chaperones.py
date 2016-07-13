import flat_files
import re
from minime.core.ProcessData import SubreactionData, PostTranslationData
from minime.core.Components import ProcessedProtein, TranslatedGene
from minime.core.MEReactions import PostTranslationReaction

# Temperature dependent Keffs for each folding mechanism/step
folding_keffs = {
    'folding_KJE_1': {'24': 0.0288, '25': 0.0323, '26': 0.0362, '27': 0.0405,
                      '28': 0.0453, '29': 0.0506, '30': 0.0565, '31': 0.0631,
                      '32': 0.0703, '33': 0.0784, '34': 0.0873, '35': 0.0971,
                      '36': 0.1080, '37': 0.1200, '38': 0.1333, '39': 0.1479,
                      '40': 0.1640, '41': 0.1818, '42': 0.2014, '43': 0.2229,
                      '44': 0.2466, '45': 0.2726, '46': 0.3012, '47': 0.3326,
                      '48': 0.3670, '49': 0.4047, '50': 0.4461},
    'folding_KJE_2': {'24': 0.0092, '25': 0.0103, '26': 0.0116, '27': 0.0130,
                      '28': 0.0146, '29': 0.0164, '30': 0.0184, '31': 0.0206,
                      '32': 0.0230, '33': 0.0258, '34': 0.0288, '35': 0.0321,
                      '36': 0.0359, '37': 0.0400, '38': 0.0446, '39': 0.0496,
                      '40': 0.0553, '41': 0.0615, '42': 0.0683, '43': 0.0759,
                      '44': 0.0842, '45': 0.0934, '46': 0.1036, '47': 0.1147,
                      '48': 0.1270, '49': 0.1405, '50': 0.1554},
    'folding_GroEL_ES': {'24': 0.2636, '25': 0.2933, '26': 0.3260,
                         '27': 0.3622,
                         '28': 0.4021, '29': 0.4461, '30': 0.4946,
                         '31': 0.5480,
                         '32': 0.6068, '33': 0.6714, '34': 0.7424,
                         '35': 0.8204,
                         '36': 0.9061, '37': 1.0000, '38': 1.1030,
                         '39': 1.2158,
                         '40': 1.3394, '41': 1.4746, '42': 1.6225,
                         '43': 1.7842,
                         '44': 1.9608, '45': 2.1537, '46': 2.3642,
                         '47': 2.5937,
                         '48': 2.8439, '49': 3.1165, '50': 3.4134}}

folding_subreactions = {
    'folding_KJE_1': {'enzymes': ['DnaJ_dim_mod_4:zn2', 'DnaK_mono'],
                      'stoichiometry': {'atp_c': -1,
                                        'h2o_c': -1,
                                        'adp_c': 1,
                                        'h_c': 1,
                                        'pi_c': 1}},
    'folding_KJE_2': {'enzymes': ['GrpE_dim'], 'stoichiometry': {}},

    'folding_GroEL_ES': {'enzymes': ['transGroES_hepta',
                                     '[GroL]14['
                                     'GroS]7_cis_with_7_adp_and_7_mg2'],
                         'stoichiometry': {'atp_c': -7,
                                           'adp_c': 7,
                                           'h_c': 7,
                                           'pi_c': 7},
                         },
    'folding_spontaneous': {'enzymes': None, 'stoichiometry': {}}
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
                     "b0002", "b2913", "b1779", "b4388", "b3919", "b1264"]


def add_chaperone_subreactions(model):
    temperature = model.global_info['temperature']
    for subreaction, values in folding_subreactions.items():
        data = SubreactionData(subreaction, model)
        data.enzyme = values['enzymes']
        data.stoichiometry = values['stoichiometry']
        if subreaction != 'folding_spontaneous':
            data.keff = folding_keffs[subreaction][str(temperature)] * 3600.


def add_chaperone_network(model):
    # First remove folding subreactions from translation reactions
    for data in model.translation_data:
        for subreaction in list(data.subreactions.keys()):
            if 'folding' in subreaction:
                data.subreactions.pop(subreaction)

    dill_df = flat_files.get_dill_keq_df()
    oobatake_df = flat_files.get_oobatake_keq_df()
    k_folding_df = flat_files.get_folding_rates_df()
    aggregation_propensity_df = flat_files.get_aggregation_popensity_df()

    for protein in model.metabolites.query(re.compile('^protein_b[0-9]')):
        if not isinstance(protein, TranslatedGene):
            continue
        protein_bnum = protein.id.replace('protein_', '')
        folded_met = ProcessedProtein(protein.id + '_folded', protein.id)
        model.add_metabolites([folded_met])

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
            if folding == 'folding_KJE_2':
                continue
            folding_id = 'folding_' + protein.id + '_' + \
                         folding.replace('_1','')

            if folding == 'folding_spontaneous':

                data = PostTranslationData(folding_id, model,
                                           protein.id + '_folded', protein.id)

                data.folding_mechanism = folding
                data.subreactions[folding] = 1
                data.aggregation_propensity = propensity
                data.k_folding = k_folding
                data.keq_folding = keq_folding


                rxn = PostTranslationReaction(folding_id)
                model.add_reaction(rxn)
                rxn.posttranslation_data = data
                rxn.update()
            else:
                data = PostTranslationData(folding_id, model,
                                            protein.id + '_folded', protein.id)

                data.subreactions[folding] = 1
                if folding == 'folding_KJE_1':
                    data.subreactions['folding_KJE_2'] = 1

                data.folding_mechanism = folding
                data.aggregation_propensity = propensity
                data.k_folding = k_folding
                data.keq_folding = keq_folding

                rxn = PostTranslationReaction(folding_id)
                model.add_reaction(rxn)
                rxn.posttranslation_data = data
                rxn.update()

        #data_2 = PostTranslationData('folding_' + protein.id + '_2', me,
        #                             protein.id + '_folded', protein.id + '_partially_folded')
        #data_2.folding_mechanism = folding
        #data_2.aggregation_propensity = propensity
        #data_2.k_folding = k_folding
        #data_2.keq_folding = keq_folding

        #rxn_2 = PostTranslationReaction('folding_' + protein.id + '_2')
        #me.add_reaction(rxn_2)
        #rxn_2.posttranslation_data = data_2
        #rxn_2.update()


def change_temperature(model, temperature):
    int_temp = int(temperature)
    str_temp = str(temperature)
    model.global_info['temperature'] = int_temp
    for subreaction, keff_dict in folding_keffs.items():
        subreaction_data = model.subreaction_data.get_by_id(subreaction)
        subreaction_data.keff = keff_dict[str_temp]

    model.update()
