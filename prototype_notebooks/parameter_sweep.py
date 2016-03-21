from cloudpickle import load, dump
import json

from minime.solve.algorithms import binary_search


def get_model():
    with open("prototype_48.pickle", "rb") as infile:
        model = load(infile)
    with open("prototype_48_expressions.pickle", "rb") as infile:
        expressions = load(infile)
    return model, expressions


def vary_unmodeled_protein_fraction(fraction):
    model, expressions = get_model()
    model.unmodeled_protein_fraction = fraction
    binary_search(model, max_mu=1.5, mu_accuracy=1e-15, verbose=True,
                  compiled_expressions=expressions)
    with open("unused_fraction_%.4f_flux.json" % fraction, "wb") as outfile:
        json.dump(model.get_metabolic_flux(), outfile)
    with open("unused_fraction_%.4f_sol.pickle" % fraction, "wb") as outfile:
        dump(model.solution, outfile)

if __name__ == "__main__":
    from multiprocessing import Pool
    from numpy import linspace

    pool = Pool(processes=4)
    for i in linspace(0, 1, 20):
        pool.apply_async(vary_unmodeled_protein_fraction, (i,))
    pool.close()
    pool.join()
