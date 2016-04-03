from parameter_sweep import *
from build_multi import *
from essentiality import *

SLURM_TEMPLATE = """#!/usr/bin/zsh
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=03:00:00
"""

CMD_TEMPLATE = """python -c "from ecolime.scripts import *; %s('%s')" """


def _get_function_name(function_name):
    # the function is passed in directly
    if hasattr(function_name, "__call__") and \
            hasattr(function_name, "__name__"):
        return function_name.__name__
    elif isinstance(function_name, string_types):
        try:
            eval(function_name)
        except:
            raise ValueError("'%s' is not the name of a function" %
                             function_name)
        else:
            return function_name
    else:
        raise TypeError("'%s' not a function or string" % function_name)


def slurm_farm(function_name, values):
    function = _get_function_name(function_name)
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


def tmux_farm(function_name, values):
    raise NotImplementedError("not yet finished")
    function_name = _get_function_name(function_name)
    for v in values:
        if not isinstance(v, string_types):
            raise ValueError("value %s is not a string" % repr(v))
        job_name = "me_%s_%s" % (function_name, v)
        command = CMD_TEMPLATE % (function_name, v)
        print("tmux new-session -s %s -d '%s' " % (job_name, command))
        system("tmux new-session -s %s -d '%s' " % (job_name, command))


def run_pool(function, values, processes=4):
    if not hasattr(function, "__call__"):
        function = _get_function_name(function)
    from multiprocessing import Pool

    pool = Pool(processes=processes)
    pool.map(function, values)
