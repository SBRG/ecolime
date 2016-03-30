from cPickle import load, dump
import re

from six import iteritems
import pandas
import sympy
from sympy.logic.boolalg import to_dnf

import cobra

from minime.solve.algorithms import binary_search, solve_at_growth_rate, \
    compile_expressions
from minime.util.building import *
from minime.core.ProcessData import *


def load_full_model():
    with open("full_model_49.pickle", "rb") as infile:
        model = load(infile)
    for exchange in model.reactions.query(re.compile("^EX_")):
        if exchange.lower_bound == 0:
            exchange.lower_bound = -10
    return model


def escape_path(name):
    return name.replace(r"/", "_SLASH_")


def create_strain_model(strain_name, model_name, homologous_loci, sequences,
                        verbose=True, solve=True):
    model = load_full_model()
    model.id = model_name + "_ME"
    model.name = strain_name
    # remove transcription
    for reaction in list(model.reactions):
        if reaction.id == "transcription_RNA_dummy":
            continue
        data = getattr(reaction, "transcription_data", None)
        if data is None:
            continue
        # we don't have tRNA/rRNA handled yet
        # TODO get files from Jon which have it
        data = reaction.transcription_data
        if len(data.RNA_products) > 0 and set(data.RNA_types) != {"mRNA"}:
            if verbose:
                print("not removing %s, which has %s" %
                      (repr(reaction), str(set(data.RNA_types))))
            continue
        model.transcription_data.remove(reaction.transcription_data.id)
        model.process_data.remove(reaction.transcription_data.id)
        reaction.remove_from_model()

    for data in list(model.transcription_data):
        continue  # this whole block is unnecessary
        if len(data.RNA_products) > 0 and set(data.RNA_types) != {"mRNA"}:
            continue
        else:
            for r in data.parent_reactions:
                r.remove_from_model()
            model.transcription_data.remove(data)
            model.process_data.remove(data)

    # remove proteins with missing homologs
    missing_proteins = set(homologous_loci.index[homologous_loci.isnull()])
    for bnum in missing_proteins:
        try:
            # for less stringent, only remove the constraint, like this
            # model.metabolites.get_by_id("protein_" + bnum).remove_from_model()
            model.reactions.get_by_id("translation_" + bnum).remove_from_model()
        except KeyError:
            if verbose:
                print("gene for removal '%s' in strain %s not found" %
                      (bnum, strain_name))
        try:
            data = model.translation_data.get_by_id(bnum)
        except KeyError:
            pass
        else:
            model.process_data.remove(data)
            model.translation_data.remove(data)

    locus_to_bnum = {locus: bnum for bnum, locus
                     in homologous_loci.dropna().iteritems()}

    for locus, na in sequences.iteritems():
        bnum = locus_to_bnum.get(locus, locus)
        na = na.upper()

        # add new transcription with fake TU's
        add_transcription_reaction(model, "TU_" + locus, {bnum}, na, update=True)

        # update translation homolog
        try:
            translation = model.translation_data.get_by_id(bnum)
        except KeyError:
            add_translation_reaction(model, bnum, dna_sequence=na, update=True)
            print("adding translation data for %s" % bnum)
        else:
            translation.nucleotide_sequence = na

    model.update()

    # add in non-k12 content
    iJO1366 = cobra.io.read_sbml_model("m_models/iJO1366.xml")
    k12_reactions = set(iJO1366.reactions.list_attr("id"))

    m_model = cobra.io.read_sbml_model("m_models/%s.xml" % model_name)
    m_reactions = set(m_model.reactions.list_attr("id"))
    specific_reactions = m_reactions.difference(k12_reactions)

    for r_id in specific_reactions:
        m_reaction = m_model.reactions.get_by_id(r_id)
        # noncatalyzed reactions
        if r_id.startswith("EX_") or r_id.startswith("DM_") or \
                "s0001" in {str(i) for i in m_reaction.genes}:
            model.add_reaction(m_reaction)
            continue
        try:
            data = model.stoichiometric_data.get_by_id(r_id)
        except KeyError:
            data = StoichiometricData(r_id, model)
            data.lower_bound = m_reaction.lower_bound
            data.upper_bound = m_reaction.upper_bound
            data._stoichiometry = {str(k): v for k, v
                                   in iteritems(m_reaction.metabolites)}

        # add in orphan reactions
        if len(m_reaction.gene_reaction_rule) == 0:
            print("adding in orphan reaction '%s'" % r_id)

            add_complexes_and_rxn_data(model, data, ["CPLX_dummy"],
                                       create_new=True, update=True)
            continue
        # ok those are the easy ones. For the rest, create a complex
        sympy_rule = sympy.sympify(m_reaction.gene_reaction_rule.replace(" and ", " & ").replace (" or ", " | "))
        converted_rule = to_dnf(sympy_rule)
        # if it's an and command, it's just one complex
        gpr_complexes = [converted_rule] if converted_rule.func == sympy.And else converted_rule.args
        complexes = []
        for count, entry in enumerate(gpr_complexes):
            if entry.is_Symbol:
                composition = [entry.name]
            elif entry.is_Boolean:
                composition = [b.name for b in entry.args]
            loci = [locus_to_bnum.get(a, a) for a in composition]
            # TODO try to use an existing complex if possible
            composition = {"protein_" + i: 1 for i in loci}
            complex_data = ComplexData("complex_%s_%d" % (r_id, count), model)
            complex_data.stoichiometry = composition
            complex_data.create_complex_formation(verbose=True)
            complexes.append(complex_data)
        add_complexes_and_rxn_data(model, data, [i.complex_id for i in complexes],
                                   create_new=True, update=True)

    model.prune()
    # no idea why this is necessary - there must be a bug somewhere
    model.reactions.biomass_dilution.reaction = \
        model.reactions.biomass_dilution.reaction

    if solve:
        expr = compile_expressions(model)
        check_sol = solve_at_growth_rate(model, 0.1, compiled_expressions=expr)
        if check_sol.status == "optimal":
            if verbose:
                print("model %s grows at 0.1" % strain_name)
            binary_search(model, verbose=verbose, compiled_expressions=expr,
                          debug=True, mu_accuracy=1e-15)
        else:
            if verbose:
                print("model %s can't grow at 0.1" % strain_name)
    with open(escape_path("me_%s.pickle" % strain_name), "wb") as outfile:
        dump(model, outfile, 2)
    return model


def get_info_frame():
    # from supplement of Jon's PNAS paper
    return pandas.read_excel("sd01.xlsx", sheetname="Table 1", skiprows=2,
                             skip_footer=1, index_col="Strain")


def get_conservation_table():
    return pandas.read_csv("homology_matrix.csv", index_col="gene")


def get_genome_sequences(genome):
    filename = "seqs/%s_seqs.csv" % genome
    seqs = pandas.read_csv(filename, index_col="locus_tag")["seq"]
    # replace all invalid characters with C
    valid_bases = re.compile("^[ATGC]$")
    invalid_bases = re.compile("[^ATGC]")
    for locus, s in seqs.iteritems():
        if not valid_bases.match(s):
            for c in invalid_bases.findall(s):
                s = s.replace(c, "C")
            seqs[locus] = s
    return seqs


def run_builder(strain_name):
    info = get_info_frame().ix[strain_name]
    model_name = info["Model Name"]
    conservation_table = get_conservation_table()
    # The genome in info doesn't have the version, and the conservation table
    # does. For example, it will be 'CU651637' in the info and 'CU651637_1' in
    # the conservation table.
    genome = info["NCBI ID"]
    genome_key = genome + "_1"

    create_strain_model(strain_name, model_name,
                        conservation_table[genome_key],
                        get_genome_sequences(genome))


def get_sources(strain_key, model_info, conservation_table, all_seqs_table):
    strain_name = strain_key.rsplit("(", 1)[0]
    info = model_info.ix[strain_name]
    model_name = info["Model Name"]

    homologous_loci = conservation_table[strain_key + "_LOCI"]
    sequences = all_seqs_table[all_seqs_table.model == model_name]
    return (strain_name, model_name, homologous_loci, sequences)


if __name__ == "__main__":
    from multiprocessing import Pool
    from os.path import isfile

    #pool = Pool(processes=4)
    # files containing homology sequences
    protein_consv = pandas.read_excel("Ecoli_k-12_gene_consv.xlsx")
    #DNA_consv = pandas.read_excel("Ecoli_k-12_gene_consv_NA_seqs.xlsx")

    all_seqs = pandas.read_csv("55_strains_seqs.csv", index_col="locus")

    # to map strain name to model name
    model_info = pandas.read_excel("sd01.xlsx", sheetname="Table 1",
                                   skiprows=2, skip_footer=1,
                                   index_col="Strain")
    model_info.pop("Unnamed: 0")
    strain_to_model = {strain: model for strain, model
                       in model_info["Model Name"].iteritems()}

    n_strains = len(protein_consv.keys()) / 4
    for i in xrange(n_strains):
        strain_key = protein_consv.keys()[4 * i + 1].rsplit("_", 1)[0]
        strain_name = strain_key.rsplit("(", 1)[0]
        model_name = strain_to_model[strain_name]

        protein_seqs = protein_consv[strain_key + "_SEQ"]
        DNA_seqs = DNA_consv[strain_key + "_SEQ"]



        if isfile(escape_path("me_%s.pickle" % strain_name)):
            print("skipped " + strain_name)
        #pool.apply_async(create_strain_model,
        #                 (strain_name, protein_seqs, DNA_seqs, False))
    #pool.close()
    #pool.join()
