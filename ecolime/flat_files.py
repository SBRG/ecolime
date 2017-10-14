from __future__ import print_function, absolute_import, division

from collections import defaultdict
import json
from os.path import dirname, join, abspath
from warnings import warn

import cobra
import pandas
from six import iteritems

import cobrame
from ecolime import corrections

ecoli_files_dir = join(dirname(abspath(__file__)), 'building_data/')

del dirname, abspath


def fixpath(filename):
    return join(ecoli_files_dir, filename)


def fix_id(id_str):
    return id_str.replace("_DASH_", "__")


def get_tu_dataframe(filename):
    tu_df = pandas.read_csv(join(ecoli_files_dir, filename), delimiter="\t",
                            index_col=0)

    tu_df = corrections.correct_tu_dataframe(tu_df)

    return tu_df


def get_complex_subunit_stoichiometry(complex_stoichiometry_file,
                                      rna_components=set()):
    """Returns dictionary of complex: {stoichiometry: {bnumber: stoichiometry}}

    some entries in the file need to be renamed.
    Colton 7/8/15 made changes directly to flat file
    renames = {"MnmE_": "b3706", "MnmG_": "b3741", "YheM_": "b3344",
    "YheL_": "b3343", "YheN_": "b3345"}

    """

    complex_stoichiometry = \
        pandas.read_table(fixpath(complex_stoichiometry_file),
                          names=['Complex', 'Name', 'Stoichiometry',
                                 'Source']).set_index('Complex')

    complex_stoichiometry_dict = {}

    for key, row in complex_stoichiometry.iterrows():
        if key.startswith('#'):
            continue

        if key in complex_stoichiometry_dict.keys():
            warn('Complex (%s) in complex_stoichiometry_file twice' % key)
        else:
            complex_stoichiometry_dict[key] = {}

        for bnums in row['Stoichiometry'].split(' AND '):

            bnum, num = bnums.rstrip(')').split('(')

            stoichiometry = float(num) if not num == '' else 1.

            prefix = 'protein_' if bnum not in rna_components else 'RNA_'
            complex_stoichiometry_dict[key][prefix + bnum] = stoichiometry

    complex_stoichiometry_dict = \
        corrections.correct_complex_stoichiometry(complex_stoichiometry_dict)

    return complex_stoichiometry_dict


def get_complex_modifications(complex_modification_file, protein_complex_file):
    """

    Reads from protein_complexes.txt and protein_modification.txt


    """

    complex_mods = pandas.read_table(fixpath(complex_modification_file))
    complex_mods = complex_mods.set_index('Modified_enzyme')

    complex_set = \
        set(get_complex_subunit_stoichiometry(protein_complex_file).keys())

    # ignore complexes which are produced in the reaction matrix
    rxn_dict = get_reaction_matrix_dict('reaction_matrix.txt',
                                        complex_set=complex_set)
    ignored_complexes = set()
    for met_stoich in rxn_dict.values():
        for met, value in iteritems(met_stoich):
            if 'mod_c' not in met:
                ignored_complexes.add(met.replace('_c', ''))
            else:
                ignored_complexes.add(met)
    # don't ignore these. They are included in the reaction matrix but still
    # must be formed via a complex formation reaction
    # TODO look into this list closer
    ignored_complexes.remove('CPLX0-782_mod_2:4fe4s')

    new_mod_dict = {}
    for key, value in iteritems(complex_mods.T.to_dict()):
        if key.startswith('#') or key in ignored_complexes:
            continue
        key = key.replace('_DASH_', '__')
        new_mod_dict[key] = {}
        new_mod_dict[key]['core_enzyme'] = value['Core_enzyme']
        new_mod_dict[key]['modifications'] = {}
        for mods in value['Modifications'].split(' AND '):
            mod, num_mods = mods.rstrip(')').split('(')
            if num_mods == '':
                num_mods = 1.
            else:
                num_mods = float(num_mods)

            mod = mod.replace('_DASH_', '__')
            new_mod_dict[key]['modifications'][mod + '_c'] = -num_mods

    new_mod_dict = corrections.correct_complex_modification_dict(new_mod_dict)

    return new_mod_dict


def get_reaction_to_complex(m_model, modifications=True):
    """anything not in this dict is assumed to be an orphan"""

    rxn_to_complex_dict = defaultdict(set)

    # Load enzyme reaction association dataframe
    df = pandas.read_csv(fixpath('enzyme_reaction_association.txt'),
                         delimiter='\t', names=['Reaction', 'Complexes'])
    # Fix legacy naming
    df = df.applymap(lambda x: x.replace('DASH', ''))
    df = df.set_index('Reaction')

    df = corrections.correct_enzyme_reaction_association_frame(df)

    for reaction, complexes in df.itertuples():
        for cplx in complexes.split(' OR '):
            if modifications:
                rxn_to_complex_dict[reaction].add(cplx)
            else:
                rxn_to_complex_dict[reaction].add(cplx.split('_mod_')[0])

    for reaction in m_model.reactions:
        if "s0001" in reaction.gene_reaction_rule:
            rxn_to_complex_dict[reaction.id].add(None)

    return rxn_to_complex_dict


def get_reaction_matrix_dict(reaction_matrix_file, complex_set=set()):
    """Return dictionary representation of the metabolic reaction matrix.
    Updates metabolite id with compartment if not contained in complex_list
    """
    matrix_df = pandas.read_csv(fixpath(reaction_matrix_file), delimiter='\t',
                                names=['Reaction', 'Metabolites',
                                       'Compartment', 'Stoichiometry'])
    matrix_df.replace({'No_Compartment': 'Cytosol'}, inplace=True)

    compartments = {'Cytosol': 'c', 'Periplasm': 'p', 'Extra-organism': 'e'}
    metabolic_reaction_dict = defaultdict(dict)
    for i, row in matrix_df.iterrows():
        reaction = fix_id(row['Reaction'])
        metabolite = fix_id(row['Metabolites'])

        # erpA is annotated incorrectly
        metabolite = metabolite.replace('CPLX-7524_mod_mn2', 'CPLX0-7617')
        stoichiometry = row['Stoichiometry']
        compartment_id = '_%s' % compartments.get(row['Compartment'])
        # use compartment to append appropriate suffix
        if metabolite.split('_mod_')[0] not in complex_set:
            metabolite += compartment_id
        metabolic_reaction_dict[reaction][metabolite] = float(stoichiometry)

    metabolic_reaction_dict = \
        corrections.correct_reaction_matrix(metabolic_reaction_dict)

    return metabolic_reaction_dict


def get_reaction_info_frame(reaction_info_file):
    df = pandas.read_csv(fixpath(reaction_info_file), delimiter="\t",
                         index_col=0)
    df = corrections.correct_reaction_info_frame(df)

    return df


def remove_compartment(id_str):
    return id_str.replace('_c', '').replace('_p', '').replace('_e', '')


def process_m_model(m_model, metabolites_file, m_to_me_map_file,
                    reaction_info_file, reaction_matrix_file,
                    protein_complex_file, defer_to_rxn_matrix=set()):

    m_model = m_model.copy()

    met_info = pandas.read_csv(fixpath(metabolites_file), delimiter="\t",
                               header=None, index_col=0,
                               names=["id", "name", "formula", "compartment",
                                      "data_source"])
    met_info.rename(lambda x: x.replace('_DASH_', '__'), inplace=True)

    complex_set = \
        set(get_complex_subunit_stoichiometry(protein_complex_file).keys())

    rxn_info = get_reaction_info_frame(reaction_info_file)
    reaction_matrix_dict = get_reaction_matrix_dict(reaction_matrix_file,
                                                    complex_set=complex_set)

    m_to_me_df = pandas.read_csv(fixpath(m_to_me_map_file), index_col=0,
                                 names=['m_name', 'me_name'])

    for rxn in list(m_model.reactions):
        if rxn.id.startswith('EX_') or rxn.id.startswith('DM_'):
            continue
        if rxn.id not in reaction_matrix_dict.keys() \
                or rxn.id in defer_to_rxn_matrix:
            rxn.remove_from_model(remove_orphans=True)
    for rxn_id in reaction_matrix_dict:
        if rxn_id in m_model.reactions:
            continue
        rxn_stoichiometry = reaction_matrix_dict[rxn_id]
        for met in rxn_stoichiometry:
            try:
                met_obj = m_model.metabolites.get_by_id(met)
            except KeyError:
                met_obj = cobrame.Metabolite(str(met))
                m_model.add_metabolites([met_obj])
            met_id = remove_compartment(met_obj.id)
            if met_id in met_info.index and not met_obj.formula:
                met_obj.formula = met_info.loc[met_id, 'formula']
                met_obj.name = met_info.loc[met_id, 'name']

        rxn = cobrame.MEReaction(rxn_id)
        m_model.add_reactions([rxn])
        rxn.add_metabolites(rxn_stoichiometry)
        reversible = rxn_info.loc[rxn_id, 'is_reversible']
        rxn.lower_bound = -1000 if reversible else 0

    for met in list(m_model.metabolites):
        met_id = remove_compartment(met.id)
        if met_id not in met_info.index and met_id in m_to_me_df.index:
            met_id = m_to_me_df.loc[met.id, 'me_name']
            if met_id != '' and met_id != 'N/A':
                met.id = met_id
            else:
                met.remove_from_model()

    # Add formulas not included in metabolites.txt
    corrections.update_metabolite_formulas(m_model)

    m_model.repair()

    return m_model


def get_m_model():
    m = cobra.Model("e_coli_ME_M_portion")
    m.compartments = {"p": "Periplasm", "e": "Extra-organism", "c": "Cytosol"}
    compartment_lookup = {v: k for k, v in iteritems(m.compartments)}

    met_info = pandas.read_csv(join(ecoli_files_dir, "metabolites.txt"),
                               delimiter="\t", header=None, index_col=0,
                               names=["id", "name", "formula", "compartment",
                                      "data_source"])
    complex_set = \
        set(get_complex_subunit_stoichiometry("protein_complexes.txt").keys())

    for met_id in met_info.index:
        fixed_id = fix_id(met_id)
        for compartment in met_info.compartment[met_id].split("AND"):
            compartment = compartment.strip()
            if compartment == "No_Compartment":
                print("Assigned %s to c" % met_id)
                compartment = m.compartments["c"]
            new_met = cobra.Metabolite(
                fixed_id + "_" + compartment_lookup[compartment])
            new_met.name = met_info.name[met_id]
            new_met.formula = met_info.formula[met_id]
            m.add_metabolites(new_met)

    rxn_info = get_reaction_info_frame('reactions.txt')
    rxn_dict = get_reaction_matrix_dict('reaction_matrix.txt',
                                        complex_set=complex_set)
    for rxn_id in rxn_info.index:
        reaction = cobrame.MEReaction(rxn_id)
        reaction.name = rxn_info.description[rxn_id]
        for met_id, amount in iteritems(rxn_dict[rxn_id]):
            try:
                metabolite = m.metabolites.get_by_id(met_id)
            except KeyError:
                metabolite = cobra.Metabolite(met_id)
            reaction.add_metabolites({metabolite: amount})
        reaction.lower_bound = \
            -1000. if rxn_info.is_reversible[rxn_id] else 0.
        reaction.upper_bound = 1000.
        if rxn_info.is_spontaneous[rxn_id]:
            reaction.gene_reaction_rule = "s0001"
        m.add_reaction(reaction)

    sources_sinks = pandas.read_csv(
        fixpath("reaction_matrix_sources_and_sinks.txt"),
        delimiter="\t", header=None, names=["rxn_id", "met_id", "compartment",
                                            "stoic"], index_col=1)

    source_amounts = pandas.read_csv(join(ecoli_files_dir,
                                          "exchange_bounds.txt"),
                                     delimiter="\t", index_col=0,
                                     names=["met_id", "amount"])

    sources_sinks.index = [fix_id(i) for i in sources_sinks.index]
    source_amounts.index = [fix_id(i) for i in source_amounts.index]

    for met in sources_sinks.index:
        met_id = met + "_" + compartment_lookup[sources_sinks.compartment[met]]
        # EX_ or DM_ + met_id
        reaction_id = sources_sinks.rxn_id[met][:3] + met_id
        reaction = cobrame.MEReaction(reaction_id)
        m.add_reaction(reaction)
        reaction.add_metabolites({m.metabolites.get_by_id(met_id): -1})
        # set bounds on exchanges
        if reaction.id.startswith("EX_") and met in source_amounts.index:
            reaction.lower_bound = -source_amounts.amount[met]

    # Add formulas not included in metabolites.txt
    corrections.update_metabolite_formulas(m)

    return m


def get_trna_modification_targets():
    trna_mod_dict = defaultdict(dict)
    filename = fixpath('post_transcriptional_modification_of_tRNA.txt')
    trna_mod = pandas.read_csv(filename, delimiter='\t')
    for mod in trna_mod.iterrows():
        mod = mod[1]
        mod_loc = '%s_at_%s' % (mod['modification'], mod['position'])
        trna_mod_dict[mod['bnum']][mod_loc] = 1

    return trna_mod_dict


def get_reaction_keffs(me, verbose=True):
    def log(*args, **kwargs):
        if verbose:
            print(*args, **kwargs)
    with open(fixpath('keffs.json'), 'r') as infile:
        keffs = json.load(infile)
    new_keffs = {}
    for r in me.reactions:
        # skip spontaneous reactions
        if getattr(r, "complex_data", None) is None:
            continue
        if isinstance(r, cobrame.MetabolicReaction) and \
                r.complex_data.id != "CPLX_dummy":
            met_rxn = r
            key = met_rxn.id.replace("-", "_DASH_").replace(
                "__", "_DASH_").replace(":", "_COLON_")

            # specific patches for PGK, TPI ids
            key = key.replace('TPI_DASH_CPLX', 'TPI')
            key = key.replace('PGK_DASH_CPLX', 'PGK')
            # key = met_rxn.id
            key = "keff_" + key.replace("_FWD_", "_").replace("_REV_", "_")

            matches = [i for i in keffs if key in i]

            # get the direction
            if met_rxn.reverse:
                matches = [i for i in matches
                           if i.endswith("_reverse_priming_keff")]
            else:
                matches = [i for i in matches
                           if i.endswith("_forward_priming_keff")]

            if len(matches) == 1:
                new_keffs[met_rxn.id] = keffs[matches[0]]
            elif len(matches) > 0:
                if len(matches) == len([i for i in matches if key + "_mod_"]):
                    new_keffs[met_rxn.id] = keffs[matches[0]]
                else:
                    log(key, len(matches))
            else:  # len(matches) == 0
                log("no keff found for " + key)
    return new_keffs


def get_m_to_me_metabolite_mapping():
    """returns a mapping from m metabolites to me metabolites"""
    f = pandas.read_csv(fixpath("m_to_me_mets.csv"), index_col=0)["me_name"]
    return f.dropna().to_dict()


def get_replace_function(source):
    def fix_columns(x):
        return x.replace(source, '').replace('C', '')
    return fix_columns


def get_dill_keq_df():
    """returns the dill length-based approximation of protein folding keqs"""
    df = pandas.read_csv(fixpath('Dill_dG_matrix.csv'))

    dill = df.rename(columns=get_replace_function('Dill_Keq_'))
    return dill.set_index('genome_region')


def get_oobatake_keq_df():
    """returns the Oobatake prediction protein folding keqs"""
    df = pandas.read_csv(fixpath('Oobatake_Keq_matrix.csv'))

    oobatake = df.rename(columns=get_replace_function('Oobatake_Keq_'))
    return oobatake.set_index('genome_region')


def get_folding_rates_df():
    """returns the Oobatake prediction protein folding keqs"""
    df = pandas.read_csv(fixpath('Folding_Rates_matrix_slope22000.csv'))

    folding_rates = df.rename(columns=get_replace_function('k_f_'))
    return folding_rates.set_index('genome_region')


def get_aggregation_popensity_df():
    """returns the Oobatake prediction protein folding keqs"""
    df = pandas.read_csv(fixpath('DnaK_reactions_parameters_5.csv'))

    return df.set_index('gene')
