import re
from collections import defaultdict

from ecolime.ecoli_k12 import *
from ecolime import ecoli_k12
import cobra
import pandas
from six import iteritems

from os.path import dirname, join, abspath

ecoli_files_dir = dirname(abspath(__file__))

del dirname, abspath


def fixpath(filename):
    return join(ecoli_files_dir, filename)


def get_complex_to_bnum_dict(rna_components):
    """Returns dictions of complex: bnumber stoichiometry

    Reads from protein_complexes.txt
    """
    ME_complex = open(fixpath('protein_complexes.txt'))
    ME_complex_dict = {}

    for line in ME_complex:
        if line.startswith("#"):
            continue
        line = line.rstrip('\tM_protein_recon\n')
        line = line.rstrip('\t2011_Updated_E_recon\n')
        line = re.split('\t| AND |', line)
        complex_composition = {}
        for gene in line[2:]:
            stoichiometry = gene[6]
            bnum = gene[0:5]
            comp_id = "RNA_" + bnum if bnum in rna_components \
                else "protein_" + bnum
            try:
                complex_composition[comp_id] = float(stoichiometry)
            except:
                complex_composition[comp_id] = float(1)
        ME_complex_dict[line[0]] = complex_composition
    ME_complex.close()
    # add in missing entries
    ME_complex_dict["CPLX0-7617"] = {"protein_b0156": 2}
    return ME_complex_dict


def get_reaction_to_complex(modifications=True, generic=False):
    """anything not in this dict is assumed to be an orphan"""
    enzRxn = open(fixpath('enzyme_reaction_association.txt'), 'r')
    rxnToModCplxDict = {}
    for line in enzRxn:
        line = line.rstrip('\n')
        line = re.split('\t| OR ', line)
        # Colton Update
        # fix legacy naming. TODO fix lysing modification with _DASH_
        if 'DASH' in line[0]:
            line[0] = line[0].replace('DASH', '')
            print 'Fixed _DASH: ', line[0]

        if generic:
            for i, cplx in enumerate(line[1:]):
                for div in ecoli_k12.divalent_list:
                    if div in cplx:
                        cplx = cplx.replace(div, 'generic_divalent')
                for mono in ecoli_k12.monovalent_list:
                    if mono in cplx:
                        cplx = cplx.replace(mono, 'generic_monovalent')
                line[i+1] = cplx

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


def get_protein_modification_dict(metabolite_list, generic=False):
    """ Get dictionary of protein modifications to components for complexes
    which don't act as metabolites in metabolic reaction.

    Return: Modified_complex: [core_complex, modification_dict]

    """
    enzMod = pandas.read_csv(fixpath('protein_modification.txt'),
                             delimiter='\t', index_col=0)
    modification_dict = {}
    for mod_complex, Modenz in enzMod.iterrows():
        if mod_complex.startswith("#"):
            continue  # commented out line
        mod_dict = defaultdict(float)

        core_complex = Modenz.Core_enzyme
        mods = Modenz.Modifications.split(' AND ')
        if mod_complex not in metabolite_list or \
                mod_complex == 'CPLX-7524_mod_mn2':
            for term in mods:
                term = term.rstrip(')')
                term = term.split('(')
                mod = term[0]
                stoich = float(term[1]) if len(term[1]) > 0 else float(1)
                if generic:
                    if mod in divalent_list:
                        mod_complex = mod_complex.replace(mod, 'generic_divalent')
                        mod = 'generic_divalent'
                    if mod in monovalent_list:
                        mod_complex = mod_complex.replace(mod, 'generic_monovalent')
                        mod = 'generic_monovalent'
                mod = mod.replace('DASH', '') + '_c'
                mod_dict[mod] += -stoich
            modification_dict[mod_complex] = [core_complex, mod_dict]

    # specific patches
    modification_dict.pop('CPLX0-246_CPLX0-1342_mod_1:SH')
    modification_dict["CPLX0-246_CPLX0-1342_mod_pydx5p"] = \
        ["CPLX0-246_CPLX0-1342", {"pydx5p_c": -1}]
    modification_dict["IscS_mod_2:pydx5p"] = \
        ["IscS", {"pydx5p_c": -2}]

    return modification_dict


def fix_id(id_str):
    return id_str.replace("_DASH_", "__")


def get_generics(generic_ions=True):
    generics = {}
    if generic_ions:
        generics["generic_divalent_c"] = [i + "_c" for i in divalent_list]
        generics["generic_monovalent_c"] = [i + "_c" for i in monovalent_list]
    return generics


def get_m_model(generic_ions=True):
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
                print "Assigned %s to c" % met_id
                compartment = m.compartments["c"]
            new_met = cobra.Metabolite(
                fixed_id + "_" + compartment_lookup[compartment])
            new_met.name = met_info.name[met_id]
            new_met.formula = met_info.formula[met_id]
            m.add_metabolites(new_met)

    if generic_ions:
        generic_ions = {"divalent": ecoli_k12.divalent_list,
                        "monovalent": ecoli_k12.monovalent_list}
    else:
        generic_ions = {}

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
        join(ecoli_files_dir, "reaction_matrix_sources_and_sinks.txt"),
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
