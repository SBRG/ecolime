import re
from collections import defaultdict

from ecolime.ecoli_k12 import *
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
    return ME_complex_dict


def get_reaction_to_modified_complex(generic=False):
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
                for div in divalent_list:
                    if div in cplx:
                        cplx = cplx.replace(div, 'generic_divalent')
                for mono in monovalent_list:
                    if mono in cplx:
                        cplx = cplx.replace(mono, 'generic_monovalent')
                line[i+1] = cplx
        rxnToModCplxDict[line[0]] = set(line[1:])
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


def get_protein_modification_dict(filename, metabolite_list, generic=False):
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

    return modification_dict

