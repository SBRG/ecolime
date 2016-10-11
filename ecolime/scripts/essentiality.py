import re
from cPickle import load
import json

from cobra.flux_analysis.deletion_worker import compute_fba_deletion
from cobrame.core.MEReactions import *
from cobrame.solve import *


def get_rich_media_model():
    with open("prototype_51.pickle", "rb") as infile:
        me = load(infile)
    for exchange in me.reactions.query(re.compile("^EX_")):
        if exchange.lower_bound == 0:
            exchange.lower_bound = -10
    return me


def get_model():
    with open("prototype_64.pickle", "rb") as infile:
        me = load(infile)
    return me


def compute_gene_essentiality(str_range):
    """str_range is "0-10" for example"""
    a_str, b_str = str_range.split("-")
    a = int(a_str)
    b = int(b_str)

    me = get_model()
    lp, solver = create_lP_at_growth_rate(me, 0.1)
    solver.solve_problem(lp)

    all_genes = me.metabolites.query(re.compile("^RNA_b[0-9]"))
    b = min(b, len(all_genes) - 1)
    target_genes = all_genes[a:b]
    results = {i.id: -1 for i in target_genes}

    def update_results(iterlimit):
        print("gene\tf\t\titerations (max %d)" % iterlimit)
        for gene_RNA_id, value in list(results.items()):
            if value > -1e-9:
                continue
            gene_RNA = me.metabolites.get_by_id(gene_RNA_id)
            # indexes should exclude transcription reactions because of TU's
            indexes = [me.reactions.index(i) for i in gene_RNA.reactions
                       if not i.id.startswith("DM") and
                       not isinstance(i, TranscriptionReaction)]
            result = compute_fba_deletion(lp, solver, me, indexes,
                                          iterlimit=iterlimit)
            if not isinstance(result, float):  # catch exceptions
                result = -1
            print("%s\t%4f\t%d" % (gene_RNA.id.split("_")[1], result,
                                   lp.numIterations))
            results[gene_RNA.id] = result
        with open("essentiality_%s.json" % str_range, "wb") as outfile:
            json.dump(results, outfile, indent=True)

    update_results(3000)
    update_results(10000)
    update_results(100000)
