import json
from os import system

from cloudpickle import load, dump
from six import string_types

from minime.solve.algorithms import binary_search

SLURM_TEMPLATE = """#!/usr/bin/zsh
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00
"""

CMD_TEMPLATE = """python -c "from parameter_sweep import *; %s('%s')" """


def get_model():
    with open("prototype_48.pickle", "rb") as infile:
        model = load(infile)
    with open("prototype_48_expressions.pickle", "rb") as infile:
        expressions = load(infile)
    return model, expressions


def save_solution(model, filename_base):
    with open(filename_base + "_flux.json", "wb") as outfile:
        json.dump(model.get_metabolic_flux(), outfile)
    with open(filename_base + "_sol.pickle", "wb") as outfile:
        dump(model.solution, outfile)


def unmodeled_protein_fraction(fraction):
    str_fraction = fraction
    fraction = float(fraction)
    model, expressions = get_model()
    model.unmodeled_protein_fraction = fraction
    binary_search(model, max_mu=1.5, mu_accuracy=1e-15, verbose=True,
                  compiled_expressions=expressions)
    save_solution(model, "unmodeled_protein_fraction_" + str_fraction)


def slurm_farm(function_name, values):
    # the function is passed in directly
    if hasattr(function_name, "__call__") and \
            hasattr(function_name, "__name__"):
        function_name = function_name.__name__
    elif isinstance(function_name, string_types):
        try:
            eval(function_name)
        except:
            raise ValueError("'%s' is not the name of a function" %
                             function_name)
    else:
        raise TypeError("'%s' not a function or string" % function_name)
    for v in values:
        if not isinstance(v, string_types):
            raise ValueError("value %s is not a string" % repr(v))
        job_name = "me_%s_%s" % (function_name, v)
        job_file = job_name + ".sl"
        with open(job_file, "w") as outfile:
            outfile.write(SLURM_TEMPLATE)
            outfile.write("#SBATCH --output=slurmout_%s\n" % job_name)
            outfile.write("#SBATCH --job-name=%s\n\n" % job_name)
            outfile.write(CMD_TEMPLATE % (function_name, v))
            outfile.write("\n")
        system("sbatch " + job_file)


if __name__ == "__main__":
    unmodeled_protein_fraction(0.1)
