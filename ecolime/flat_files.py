from Bio import SeqIO
import re


divalent_list = ['ca2', 'mg2', 'mn2', 'cobalt2', 'ni2', 'cd2', 'zn2']
monovalent_list = ['k', 'na1']


# TODO: Popuate transcription reactions using function here
def get_K12_genbank():
    gb_file = SeqIO.read('NC_000913.2.gb', 'gb')
    return gb_file


def get_complex_to_bnum_dict():
    """Returns dictions of complex: bnumber stoichiometry

    Reads from protein_complexes.txt
    """
    ME_complex = open('protein_complexes.txt')
    ME_complex_dict = {}

    for line in ME_complex:
        line = line.rstrip('\tM_protein_recon\n')
        line = line.rstrip('\t2011_Updated_E_recon\n')
        line = re.split('\t| AND |', line)
        ME_complex_dict[line[0]] = line[2:]
    return ME_complex_dict
    ME_complex.close()


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
                        line[i+1] = cplx.replace(div, 'generic_divalent')
                for mono in monovalent_list:
                    if mono in cplx:
                        line[i+1] = cplx.replace(mono, 'generic_monovalent')
        rxnToModCplxDict[line[0]] = line[1:]
    return rxnToModCplxDict


def get_reaction_matrix_dict():
    reaction_matrix = open('reaction_matrix.txt', 'r')
    ME_reaction_dict = {}
    for line in reaction_matrix:
        line = line.replace('\n', '')
        line = line.split('\t')
        line[1] = line[1].replace('DASH', '')
        # use compartment to append appropriate suffix
        if line[2] == 'Cytosol':
            line[1] += '_c'
        elif line[2] == 'Periplasm':
            line[1] += '_p'
        elif line[2] == 'Extra-organism':
            line[1] += '_e'
        if line[0] not in ME_reaction_dict:
            ME_reaction_dict[line[0]] = []
        ME_reaction_dict[line[0]] += [[line[1], float(line[3])]]
        return ME_reaction_dict


def get_reaction_info_dict():
    """ Index [1] in reaction_info_dict list refers to reversibility 1 = reversible

    """

    ME_reactions = open('reactions.txt', 'r')
    Reaction_info_dict = {}
    for line in ME_reactions:
        line = line.split('\t')
        Reaction_info_dict[line[0]] = [line[1], line[2]]

    return Reaction_info_dict


def get_iJO_update_list():
    ME_reactions = open('reactions.txt', 'r')
    iJO_update_list = []
    for line in ME_reactions:
        line = line.split('\t')
        if line[3] != 'iJO1366':
            iJO_update_list.append(line[0])
    iJO_update_list = iJO_update_list[1:]
    return iJO_update_list
