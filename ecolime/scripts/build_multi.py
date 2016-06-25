from cPickle import load, dump
import re
import json

from six import iteritems
import pandas
import sympy
from sympy.logic.boolalg import to_dnf
import numpy as np

import cobra

from minime.solve.algorithms import binary_search, solve_at_growth_rate, \
    compile_expressions
from minime.util.building import *
from minime.core.ProcessData import *
from ecolime.flat_files import get_m_to_me_metabolite_mapping


def load_full_model():
    with open("full_model_54.pickle", "rb") as infile:
        model = load(infile)
    for exchange in model.reactions.query(re.compile("^EX_")):
        if exchange.lower_bound == 0:
            exchange.lower_bound = -10
    for m_id, me_id in get_m_to_me_metabolite_mapping().items():
        mapping = cobra.Reaction("rename_%s_to_%s" % (m_id, me_id))
        model.add_reaction(mapping)
        mapping.add_metabolites({
            cobra.Metabolite(m_id): -1,
            me_id: 1})
        mapping.lower_bound = -1000

    # New reaction for thf_c auxotrophs from Monk et al
    print 'adding thf_c exchange'
    rxn = cobra.Reaction('EX_thf_c')
    model.add_reaction(rxn)
    rxn.add_metabolites({model.metabolites.thf_c: -1})
    rxn.upper_bound = 1000
    rxn.lower_bound = -1000

    return model

locus_prefix_to_model_id = {'EcDH1': 'iEcDH1_1363'}


def add_homologous_protein(model, locus, info_frame):
    model_id = locus_prefix_to_model_id.get(locus.split('_')[0])
    if not model_id:
        return False
    NCBI_id = info_frame.ix[model_id]['NCBI ID']
    seqs = get_genome_sequences(NCBI_id)
    try:
        seq = seqs[locus]
    except KeyError:
        return False

    if 'RNA_' + locus not in model.metabolites:
        transcript = TranscribedGene('RNA_' + locus)
        # print("adding transcript (%s) and demand to model" % transcript.id)
        transcript.nucleotide_sequence = seq
        model.add_metabolites([transcript])
        demand_reaction = Reaction("DM_" + transcript.id)
        model.add_reaction(demand_reaction)
        demand_reaction.add_metabolites({transcript.id: -1})

    if 'transcription_TU_' + locus not in model.reactions:
        # add new transcription with fake TU's
        add_transcription_reaction(model, "TU_" + locus, {locus},
                                   seq, update=True)
        model.transcription_data.get_by_id("TU_" + locus).RNA_polymerase = \
            'RNAP70-CPLX'

    # update translation homolog
    try:
        translation = model.translation_data.get_by_id(locus)
    except KeyError:
        add_translation_reaction(model, locus, dna_sequence=seq, update=True)
        print("adding translation data for %s" % locus)
    else:
        translation.nucleotide_sequence = seq

    return True


def create_strain_model(strain_name, model_name, homologous_loci, sequences,
                        verbose=True, solve=True):

    # get essential
    with open("essentiality.json", "rb") as infile:
        essential = set(json.load(infile))

    # remove proteins with missing homologs
    missing_proteins = set(homologous_loci.index[homologous_loci.isnull()])
    missing_essential = missing_proteins.intersection(essential)
    missing_proteins.difference_update(essential)  # don't remove essential
    missing_essential_renamed = {"RNA_" + i for i in missing_essential}
    with open("%s_missing.json" % model_name, "wb") as outfile:
        json.dump(sorted(missing_essential), outfile, indent=True)

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
        shouldnt_remove = data.RNA_products.intersection(
            missing_essential_renamed)
        if len(shouldnt_remove) > 0:
            if verbose:
                print("not removing %s, which has essential genes %s" %
                      (repr(reaction), str(shouldnt_remove)))
            continue
        model.transcription_data.remove(reaction.transcription_data.id)
        model.process_data.remove(reaction.transcription_data.id)
        reaction.remove_from_model()

#    for data in list(model.transcription_data):
#        continue  # this whole block is unnecessary
#        if len(data.RNA_products) > 0 and set(data.RNA_types) != {"mRNA"}:
#            continue
#        else:
#            for r in data.parent_reactions:
#                r.remove_from_model()
#            model.transcription_data.remove(data)
#            model.process_data.remove(data)

    # remove proteins with missing homologs
    for bnum in missing_proteins:
        try:
            # for less stringent, only remove the constraint, like this
            # model.metabolites.get_by_id("protein_" + bnum).remove_from_model()
            model.reactions.get_by_id("translation_" + bnum).remove_from_model()
            for rxn in model.reactions.query('translocation_' + bnum):
                print("removing %s" % rxn.id)
                rxn.remove_from_model()
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
            for post_trans_data in model.posttranslation_data.query(bnum):
                model.posttranslation_data.remove(post_trans_data.id)

    locus_to_bnum = {locus: bnum for bnum, locus
                     in homologous_loci.dropna().iteritems()}

    for locus, na in sequences.iteritems():
        if type(na) != str:
            from IPython import embed
            embed()
        bnum = locus_to_bnum.get(locus, locus)
        na = na.upper()

        # Create new transcript and demand reaction (Must be done in
        # multi builder script _Sequence required to obtain mass)
        # Colton addition

        if 'RNA_' + bnum not in model.metabolites:
            transcript = TranscribedGene('RNA_' + bnum)
            #print("adding transcript (%s) and demand to model" % transcript.id)
            transcript.nucleotide_sequence = na
            model.add_metabolites([transcript])
            demand_reaction = Reaction("DM_" + transcript.id)
            model.add_reaction(demand_reaction)
            demand_reaction.add_metabolites({transcript.id: -1})

        if 'transcription_TU_' + locus not in model.reactions:
            # add new transcription with fake TU's
            add_transcription_reaction(model, "TU_" + locus, {bnum},
                                       na, update=True)
            model.transcription_data.get_by_id("TU_" + locus).RNA_polymerase = \
                'RNAP70-CPLX'

        # update translation homolog
        try:
            translation = model.translation_data.get_by_id(bnum)
        except KeyError:
            add_translation_reaction(model, bnum, dna_sequence=na, update=True)
            print("adding translation data for %s" % bnum)
        else:
            translation.nucleotide_sequence = na

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
            if data.lower_bound < 0:
                add_metabolic_reaction_to_model(model, data.id, 'reverse',
                                                "CPLX_dummy", update=True)
            if data.upper_bound > 0:
                add_metabolic_reaction_to_model(model, data.id, 'forward',
                                                "CPLX_dummy", update=True)
            continue
        # ok those are the easy ones. For the rest, create a complex
        sympy_rule = sympy.sympify(m_reaction.gene_reaction_rule.replace(" and ", " & ").replace (" or ", " | "))
        converted_rule = to_dnf(sympy_rule)
        # just one gene
        if isinstance(converted_rule, sympy.Symbol):
            gpr_complexes = [converted_rule]
        # if it's an and command, it's just one complex
        elif converted_rule.func == sympy.And:
            gpr_complexes = [converted_rule]
        else:  # or
            gpr_complexes = converted_rule.args

        complexes = []
        info_frame = get_info_frame()
        for count, entry in enumerate(gpr_complexes):
            if entry.is_Symbol:
                composition = [entry.name]
            elif entry.is_Boolean:
                composition = [b.name for b in entry.args]
            # TODO try to use an existing complex if possible
            fixed_loci = []
            for locus in [locus_to_bnum.get(a, a) for a in composition]:
                # Specific patched for iECDH1ME8569 TODO move these somewhere else
                if 'protein_' + locus not in model.metabolites:
                    added = add_homologous_protein(model, locus, info_frame)
                    if not added:
                        print("No %s homolog found, using dummy" % locus)
                        locus = 'dummy'
                fixed_loci.append(locus)
            if len(fixed_loci) == 0:
                raise('Complex %s must contain a protein' % "complex_" + r_id)
            composition = {"protein_" + i: 1 for i in fixed_loci}
            complex_data = ComplexData("complex_%s_%d" % (r_id, count), model)
            complex_data.stoichiometry = composition
            complex_data.create_complex_formation(verbose=True)
            complexes.append(complex_data)
        for cplx in complexes:
            if data.lower_bound < 0:
                add_metabolic_reaction_to_model(model, data.id, 'reverse',
                                                cplx.id, update=True)
            if data.upper_bound > 0:
                add_metabolic_reaction_to_model(model, data.id, 'forward',
                                                cplx.id, update=True)
    model.update()
    model.prune()

    # ultra easy mode - remove biomass components
    # for i in list(model.reactions.biomass_dilution.metabolites):
    #    if i.id != "biomass":
    #        model.reactions.biomass_dilution.pop(i)

    if solve:
        expr = compile_expressions(model)
        check_sol = solve_at_growth_rate(model, 0.1, compiled_expressions=expr)
        if check_sol.status == "optimal":
            if verbose:
                print("model %s grows at 0.1" % strain_name)
            check_0 = False
            binary_search(model, verbose=verbose, compiled_expressions=expr,
                          debug=True, mu_accuracy=1e-6)
        else:
            if verbose:
                check_0 = True
                print("model %s can't grow at 0.1" % strain_name)

        if check_0:
            check_0 = solve_at_growth_rate(model, 0, compiled_expressions=expr)
            if check_0.status == "optimal":
                if verbose:
                    print("model %s grows at 0" % strain_name)
                check_m_biomass = False
            else:
                if verbose:
                    print("model %s can't grow at 0" % strain_name)
                check_m_biomass = True

    with open("me_%s.pickle" % model_name, "wb") as outfile:
        dump(model, outfile, 2)
    return model


def get_info_frame():
    # from supplement of Jon's PNAS paper
    return pandas.read_excel("sd01.xlsx", sheetname="Table 1", skiprows=2,
                             skip_footer=1, index_col="Model Name")


def get_conservation_table():
    return pandas.read_csv("homology_matrix.csv", index_col="gene")


def get_genome_sequences(genome):
    filename = "seqs/%s_seqs.csv" % genome
    fixed_seqs = {}
    seqs = pandas.read_csv(filename).set_index("locus_tag")["seq"]
    # replace all invalid characters with C
    valid_bases = re.compile("^[ATGC]$")
    invalid_bases = re.compile("[^ATGC]")
    for locus, s in seqs.iteritems():
        if type(s) == float:
            print locus, 'has no sequence'
        elif type(locus) == float:
            print s, 'has no locus id'
        else:
            if not valid_bases.match(s):
                for c in invalid_bases.findall(s):
                    s = s.replace(c, "C")
            fixed_seqs[locus] = s

    return fixed_seqs


def run_builder(model_name):
    info = get_info_frame().ix[model_name]
    strain_name = info["Strain"]
    conservation_table = get_conservation_table()
    # The genome in info doesn't have the version, and the conservation table
    # does. For example, it will be 'CU651637' in the info and 'CU651637_1' in
    # the conservation table.
    genome = info["NCBI ID"]
    genome_key = genome + "_1"

    create_strain_model(strain_name, model_name,
                        conservation_table[genome_key],
                        get_genome_sequences(genome))
