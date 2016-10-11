from __future__ import print_function
import re
from collections import defaultdict
import json
from os.path import dirname, join, abspath


from cobrame.core.MEReactions import MetabolicReaction
import cobra
import pandas
from six import iteritems

from ecolime.ecoli_k12 import *
from ecolime import ecoli_k12, tRNA_charging

ecoli_files_dir = dirname(abspath(__file__))

del dirname, abspath


def fixpath(filename):
    return join(ecoli_files_dir, filename)


def get_complex_composition(rna_components):
    """Returns dictionary of complex: {stoichiometry: {bnumber: stoichiometry},
                                       modifications: {modificaiton: number}

    Reads from protein_complexes.txt and protein_modification.text

    some entries in the file need to be renamed.
    Colton 7/8/15 made changes directly to flat file
    renames = {"MnmE_": "b3706", "MnmG_": "b3741", "YheM_": "b3344",
    "YheL_": "b3343", "YheN_": "b3345"}
    """
    complex_stoich = pandas.read_table(fixpath('protein_complexes.txt'),
                                       names=['Complex', 'Name',
                                              'Stoichiometry',
                                              'Source']).set_index('Complex')

    complex_mods = pandas.read_table(fixpath('protein_modification.txt'))
    complex_mods = complex_mods.set_index('Modified_enzyme')

    # ignore complexes which are produced in the reaction matrix
    rxn_dict = get_reaction_matrix_dict()
    ignored_complexes = set()
    for met_stoich in rxn_dict.values():
        for met in met_stoich:
            if 'mod_c' not in met:
                ignored_complexes.add(met.replace('_c',''))
            else:
                ignored_complexes.add(met)
    # don't ignore these. They are included in the reaction matrix but still
    # must be formed via a complex formation reaction
    # TODO look into this list closer
    ignored_complexes.remove('EG10544-MONOMER_mod_palmitate')
    ignored_complexes.remove('CPLX-7524_mod_mn2')
    ignored_complexes.remove('CPLX0-782_mod_2:4fe4s')
    ignored_complexes.remove('DSBC-CPLX_mod_Oxidized')
    ignored_complexes.remove('DSBG-CPLX_mod_Oxidized')
    ignored_complexes.remove('EG11597-MONOMER_mod_amp')

    new_mod_dict = {}
    for key, value in complex_mods.T.to_dict().items():
        if key.startswith('#') or key in ignored_complexes:
            continue
        new_mod_dict[key] = {}
        new_mod_dict[key]['core_enzyme'] = value['Core_enzyme']
        new_mod_dict[key]['modifications'] = {}
        for mods in value['Modifications'].split(' AND '):
            mod, num_mods = mods.rstrip(')').split('(')
            if num_mods == '':
                num_mods = 1.
            else:
                num_mods = float(num_mods)

            mod = mod.replace('DASH', '')
            new_mod_dict[key]['modifications'][mod + '_c'] = -num_mods

    # specific patches. Not used in iOL1650 but included as a complex
    for i in ['CPLX0-246_CPLX0-1342_mod_1:SH']:
        new_mod_dict.pop(i)

    new_mod_dict["CPLX0-246_CPLX0-1342_mod_pydx5p"] = {"core_enzyme":
                                                           "CPLX0-246_CPLX0-1342",
                                                       "modifications": {"pydx5p_c": 1}}
    new_mod_dict["IscS_mod_2:pydx5p"] = {"core_enzyme": "IscS",
                                         "modifications": {"pydx5p_c": 2}}


    new_stoich_dict = {}
    for key, value in complex_stoich['Stoichiometry'].T.to_dict().items():
        if key.startswith('#'):
            continue
        new_stoich_dict[key] = {}
        for bnums in value.split(' AND '):

            bnum, num = bnums.rstrip(')').split('(')

            if num == '':
                num = 1.
            else:
                num = float(num)

            prefix = 'protein_' if bnum not in rna_components else 'RNA_'
            new_stoich_dict[key][prefix + bnum] = num

    # add in missing entries
    new_stoich_dict["CPLX0-7617"] = {"protein_b0156": 2}
    # should actually be a dimer PMID 24914049
    new_stoich_dict['YdaO_mono'] = {"protein_b1344": 2}

    return new_stoich_dict, new_mod_dict


def get_reaction_to_complex(modifications=True):
    """anything not in this dict is assumed to be an orphan"""
    enzRxn = open(fixpath('enzyme_reaction_association.txt'), 'r')
    rxnToModCplxDict = {}
    for line in enzRxn:
        line = line.rstrip('\n')
        line = re.split('\t| OR ', line)
        # Colton Update
        # fix legacy naming. TODO fix lysine modification with _DASH_
        if 'DASH' in line[0]:
            line[0] = line[0].replace('DASH', '')
            print('Fixed _DASH: ' + line[0])

        reaction = line[0]
        rxnToModCplxDict[reaction] = set()
        for line in line[1:]:
            if modifications:
                rxnToModCplxDict[reaction].add(line)
            else:
                rxnToModCplxDict[reaction].add(line.split('_mod_')[0])
    enzRxn.close()
    m_model = get_m_model()
    for reaction in m_model.reactions:
        if reaction.gene_reaction_rule == "s0001":
            if reaction.id not in rxnToModCplxDict:
                rxnToModCplxDict[reaction.id] = set()
            rxnToModCplxDict[reaction.id].add(None)
    return rxnToModCplxDict


def get_reaction_matrix_dict():
    reaction_matrix = open(fixpath('reaction_matrix.txt'), 'r')
    # These metabolites are mistakenly labeled as NoCompartment when they
    # should really be in the cytosol.
    move_to_cytosol = {'adp', 'atp', 'h', 'pi', '2tpr3dpcoa', 'dpm', 'fe2',
                       'dad__5', 'met__L', 'tl'}
    ME_reaction_dict = {}
    for line in reaction_matrix:
        line = line.strip()
        if line.startswith("#") or len(line) == 0:
            continue
        rxn, met, comp, count = line.split('\t')
        rxn = rxn.replace('DASH', '')
        met = met.replace('DASH', '')
        # use compartment to append appropriate suffix
        if comp == 'Cytosol':
            met += '_c'
        elif comp == 'Periplasm':
            met += '_p'
        elif comp == 'Extra-organism':
            met += '_e'
        # some mistakenly annotated as no compartment
        elif comp == 'No_Compartment' and met in move_to_cytosol:
            met += '_c'
        if rxn not in ME_reaction_dict:
            ME_reaction_dict[rxn] = {}
        ME_reaction_dict[rxn][met] = float(count)
    reaction_matrix.close()
    for rxn_id in ["PFL_act", "hemeD_synthesis", "23bpg_generation"]:
        ME_reaction_dict[rxn_id] = \
            {k + "_c": v for k, v in iteritems(ME_reaction_dict[rxn_id])}
    return ME_reaction_dict


def get_reaction_info_frame():
    return pandas.read_csv(fixpath("reactions.txt"),
                           delimiter="\t", index_col=0)


def fix_id(id_str):
    return id_str.replace("_DASH_", "__")


def get_m_model():
    m = cobra.Model("e_coli_ME_M_portion")
    m.compartments = {"p": "Periplasm", "e": "Extra-organism", "c": "Cytosol"}
    compartment_lookup = {v: k for k, v in m.compartments.items()}

    met_info = pandas.read_csv(join(ecoli_files_dir, "metabolites.txt"),
                               delimiter="\t", header=None, index_col=0,
                               names=["id", "name", "formula", "compartment",
                                      "data_source"])

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

    rxn_info = get_reaction_info_frame()
    rxn_dict = get_reaction_matrix_dict()
    for rxn_id in rxn_info.index:
        reaction = cobra.Reaction(rxn_id)
        reaction.name = rxn_info.description[rxn_id]
        for met_id, amount in rxn_dict[rxn_id].items():
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
        reaction = cobra.Reaction(reaction_id)
        m.add_reaction(reaction)
        reaction.add_metabolites({m.metabolites.get_by_id(met_id): -1})
        # set bounds on exchanges
        if reaction.id.startswith("EX_") and met in source_amounts.index:
            reaction.lower_bound = -source_amounts.amount[met]

    return m


def get_tRNA_modification_targets():
    tRNA_mod_dict = defaultdict(dict)
    filename = fixpath('post_transcriptional_modification_of_tRNA.txt')
    tRNA_mod = pandas.read_csv(filename, delimiter='\t')
    for mod in tRNA_mod.iterrows():
        mod = mod[1]
        mod_loc = '%s_at_%s' % (mod['modification'], mod['position'])
        tRNA_mod_dict[mod['bnum']][mod_loc] = 1

    return tRNA_mod_dict


def get_reaction_keffs(me, verbose=True):
    def log(*args, **kwargs):
        if verbose:
            print(*args, **kwargs)
    with open(fixpath('keffs.json'), 'rb') as infile:
        keffs = json.load(infile)
    new_keffs = {}
    for r in me.reactions:
        # skip spontaneous reactions
        if getattr(r, "complex_data", None) is None:
            continue
        if isinstance(r, MetabolicReaction) and r.complex_data.id != "CPLX_dummy":
            met_rxn = r
            key = met_rxn.id.replace("-", "_DASH_").replace("__", "_DASH_").replace(":", "_COLON_")
            # specific patches for PGK, TPI ids
            key = key.replace('TPI_DASH_CPLX', 'TPI')
            key = key.replace('PGK_DASH_CPLX', 'PGK')
            # key = met_rxn.id
            key = "keff_" + key.replace("_FWD_", "_").replace("_REV_", "_")

            matches = [i for i in keffs if key in i]
            # get the direction
            if met_rxn.reverse:
                matches = [i for i in matches if i.endswith("_reverse_priming_keff")]
            else:
                matches = [i for i in matches if i.endswith("_forward_priming_keff")]
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


def get_tRNA_modification_procedures():

    mod = tRNA_charging.trna_modification.copy()

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

    # no reference for this, but it's free anyways. doesn't make sense to have
    # it in, and it's probably not real anyways
    correct_mod = mod['mnm5s2U_at_34']['carriers']
    correct_mod.pop('trdrd_c')
    correct_mod.pop('trdox_c')

    return mod


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