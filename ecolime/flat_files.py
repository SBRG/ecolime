import re

import pandas

from os.path import dirname, join, abspath

ecoli_files_dir = dirname(abspath(__file__))

del dirname, abspath


divalent_list = ['ca2', 'mg2', 'mn2', 'cobalt2', 'ni2', 'cd2', 'zn2']
monovalent_list = ['k', 'na1']


def get_complex_to_bnum_dict():
    """Returns dictions of complex: bnumber stoichiometry

    Reads from protein_complexes.txt
    """
    ME_complex = open(join(ecoli_files_dir, 'protein_complexes.txt'))
    ME_complex_dict = {}

    for line in ME_complex:
        line = line.rstrip('\tM_protein_recon\n')
        line = line.rstrip('\t2011_Updated_E_recon\n')
        line = re.split('\t| AND |', line)
        ME_complex_dict[line[0]] = line[2:]
    ME_complex.close()
    return ME_complex_dict


def get_reaction_to_modified_complex(generic=False):
    enzRxn = open('enzyme_reaction_association.txt', 'r')
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
    reaction_matrix = open(join(ecoli_files_dir, 'reaction_matrix.txt'), 'r')
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
        if rxn not in ME_reaction_dict:
            ME_reaction_dict[rxn] = {}
        ME_reaction_dict[rxn][met] = float(count)
    reaction_matrix.close()
    return ME_reaction_dict


def get_reaction_info_frame():
    return pandas.read_csv("reactions.txt", delimiter="\t", index_col=0)


def get_protein_modification_dict(generic=False):
    """ Get dictionary of protein modifications to components for complexes
    which don't act as metabolites in metabolic reaction.

    Return: Modified_complex: [core_complex, modification_dict]

    """
    reaction_matrix = open(join(ecoli_files_dir, 'reaction_matrix.txt'), 'r')
    metlist = []
    for line in reaction_matrix:
        line = line.replace('\n', '')
        line = line.split('\t')
        line[1] = line[1].replace('DASH', '')
        metlist.append(line[1])

    enzMod = open(join(ecoli_files_dir, 'protein_modification.txt'), 'r')
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
    return modification_dict
