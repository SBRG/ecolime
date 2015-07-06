import re

import pandas
from six import iteritems

from os.path import dirname, join, abspath

ecoli_files_dir = dirname(abspath(__file__))

del dirname, abspath

def fixpath(filename):
    return join(ecoli_files_dir, filename)


def get_complex_to_bnum_dict():
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
        ME_complex_dict[line[0]] = line[2:]
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


def get_protein_modification_dict(generic=False):
    """ Get dictionary of protein modifications to components for complexes
    which don't act as metabolites in metabolic reaction.

    Return: Modified_complex: [core_complex, modification_dict]

    """
    reaction_matrix = open(fixpath('reaction_matrix.txt'), 'r')
    metlist = []
    for line in reaction_matrix:
        line = line.replace('\n', '')
        line = line.split('\t')
        line[1] = line[1].replace('DASH', '')
        metlist.append(line[1])

    enzMod = open(fixpath('protein_modification.txt'), 'r')
    modification_dict = {}
    for line in enzMod:
        if line.startswith("#"):
            continue  # commented out line
        mod_dict = {}
        line = line.rstrip('\tM_protein_recon\n')
        line = line.rstrip('\t2011_Updated_E_recon\n')
        line = re.split('\t| AND |', line)
        # CPLX-7524_mod_mn2 acts as reactant in FE-S reactions and as enyzme
        # for ALLTAMH
        if line[0] not in metlist or line[0] == 'CPLX-7524_mod_mn2':
            core_complex = line[1]
            for term in line[2:]:
                term = term.rstrip(')')
                term = term.split('(')
                try:
                    stoich = float(term[1])
                    mod = term[0]
                except:
                    stoich = 1
                    mod = term[0]
                if generic:
                    if mod in divalent_list:
                        line[0] = line[0].replace(mod, 'generic_divalent')
                        mod = 'generic_divalent'
                    if mod in monovalent_list:
                        line[0] = line[0].replace(mod, 'generic_monovalent')
                        mod = 'generic_monovalent'
                mod = mod.replace('DASH', '') + '_c'
                if mod not in mod_dict:
                    mod_dict[mod] = 0
                mod_dict[mod] += -stoich
            modification_dict[line[0]] = [core_complex, mod_dict]
    reaction_matrix.close()
    enzMod.close()
    # specific patches
    modification_dict.pop('CPLX0-246_CPLX0-1342_mod_1:SH')
    modification_dict["CPLX0-246_CPLX0-1342_mod_pydx5p"] = \
        ["CPLX0-246_CPLX0-1342", {"pydx5p_c": -1}]
    return modification_dict
