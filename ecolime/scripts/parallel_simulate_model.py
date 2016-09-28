# ==========================================================
# Sensitivity of substrate uptake hierarchy to ME params
#
# Laurence Yang, SBRG, UCSD
#
# 25 Mar 2016: first version
# 08 Apr 2016: sample starting at calibrated model (keffs)
# 17 Aug 2016: modified by Colton Lloyd

# =============== Supported simulation types =================
# carbon_substrates: all growth supporting carbon sources

# ============================================================
from mpi4py import MPI
from qminos.me1 import ME_NLP1
import numpy as np
import cPickle
import os
import json
import re
import time
import pandas as pd
import argparse

# ------------------------------------------------------------
# Modules to Manipulate Model
# ------------------------------------------------------------
from manipulate_model import (change_substrate, change_PO_ratio,
                              modify_model_parameter)

# ************************************************************
# Parameters
# ------------------------------------------------------------
# Bisection parameters
MU_PREC = 1e-3
MU_MIN = 0.
MU_MAX = 2.0

# Simulation parameters
ME_PROTO_VER = 62

simstr = ''

# ************************************************************
# Argument
parser = argparse.ArgumentParser(description='Simulation parameters.')

parser.add_argument('sim_type', help='', type=str)
parser.add_argument('--N_sims', help='', type=int, default=0)
parser.add_argument('--Start',  help='', default=0.)
parser.add_argument('--Stop', help='', default=0.)
parser.add_argument('--Mult',  help='', default=0.)
args = parser.parse_args()


SIMULATION_TYPE = args.sim_type
N_SIMULATIONS = int(args.N_sims)
START_VALUE = float(args.Start)
STOP_VALUE = float(args.Stop)
MULT_VALUE = float(args.Mult)


# ************************************************************
# Load and prepare model for simulations
with open('~/ecolime/prototype_notebooks/prototype_%s.pickle' % ME_PROTO_VER, 'r') as f:
    model = cPickle.load(f)

# ************************************************************

# Get MPI worker information
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

# number of cores running process
nWorkers = size
# ------------------------------------------------------------


def return_substrate_list(model, source='C'):
    """source must be C, P, S, N for carbon, phosphorous, etc"""
    rxns = []
    for rxn in model.reactions.query(re.compile('^EX_')):
        try:
            if source in rxn._metabolites.keys()[0].formula:
                rxns.append(rxn.id)
        except TypeError:
            print(rxn.id, ' has no formula')
    return rxns

if SIMULATION_TYPE == 'carbon_substrates':
    reaction_list = return_substrate_list(model, source='C')
    N_SIMULATIONS = len(reaction_list)
elif SIMULATION_TYPE == 'PO_Ratio':
    pass
elif SIMULATION_TYPE == 'glucose_growth_curve':
    pass
else:
    raise('Simulation type %s not supported' % SIMULATION_TYPE)

output_dir = SIMULATION_TYPE

if not os.path.isdir(os.path.join(os.getcwd(), output_dir)):
    os.mkdir(os.path.join(os.getcwd(), output_dir))


# ------------------------------------------------------------
# Dictionary of work:
samples = range(0, N_SIMULATIONS)

# Used for parameter sweeps / growth curves
value_list = np.linspace(START_VALUE, STOP_VALUE, N_SIMULATIONS)

work_dict_list = [{'sample': sample} for sample in samples]
inds_all = range(0, len(work_dict_list))
# -----------------------------------------------------------

###
print(nWorkers, 'workers available')
print('%d tasks' % len(inds_all))
###

# ------------------------------------------------------------
# Split the work
nTasks = len(inds_all)
tasks_per_worker = nTasks/nWorkers
rem = nTasks - tasks_per_worker * nWorkers
# Distributes remainder across workers
worker_sizes = np.array([tasks_per_worker + (1 if i < rem else 0)
                         for i in range(0, nWorkers)], 'i')

# ------------------------------------------------------------
# Thus, worker_tasks[rank] gives the actual indices to work on
inds_end = np.cumsum(worker_sizes)
inds_start = np.hstack((0, inds_end[0:-1]))
worker_tasks = [inds_all[inds_start[i]:inds_end[i]] for i in range(0,
                                                                   nWorkers)]

print('inds_all:', inds_all)
print('worker_tasks:', worker_tasks)
print('%d batches' % len(worker_tasks))

# ============================================================
# Work performed by each worker

# Compile Expressions
me_nlp = ME_NLP1(model, growth_key='mu')
me_nlp.compiled_expressions = me_nlp.compile_expressions()

def do_work(inds):
    # Solve
    def simulate(ind, hs):
        # ----------------------------------------------------
        # Choose sample
        sim_dict = work_dict_list[ind]
        sample = sim_dict['sample']
        print('Simulation sample %d in job %s' % (sample, rank))
        # ----------------------------------------------------
        tic = time.time()

        def solve_model(model):
            # Re-compile expressions with new keffs in S and solve
            me_nlp = ME_NLP1(model, growth_key='mu')
            return me_nlp.bisectmu(precision=MU_PREC, mumin=MU_MIN, mumax=MU_MAX,
                                   basis=hs)

        value = value_list[sample]
        if SIMULATION_TYPE == 'carbon_substrates':
            substrate_rxn = reaction_list[sample]
            change_substrate(model, substrate_rxn)
            output_file = '%s' % substrate_rxn
            print(output_file)
            muopt, hs, xopt, cache = solve_model(model)
            model.reactions.get_by_id(substrate_rxn).lower_bound = 0
        elif SIMULATION_TYPE == 'PO_Ratio':
            output_file = '%.2f_fraction_ndh1' % value
            change_PO_ratio(model, value)
            muopt, hs, xopt, cache = solve_model(model)
        elif 'growth_curve' in SIMULATION_TYPE:
            substrate_dict = {'glucose': 'EX_glc__D_e'}
            substate_name = SIMULATION_TYPE.split('_')[0]
            uptake_rxn = substrate_dict.get(substate_name, 'EX_glc__D_e')
            if 'parameter' in SIMULATION_TYPE:
                output_file = '%05.2f_mult_%06.3f' % (MULT_VALUE, value)
            else:
                output_file = '%06.3f_%s' % (value, uptake_rxn)
            change_substrate(model, uptake_rxn, value=value)
            muopt, hs, xopt, cache = solve_model(model)
        else:
            raise('Simulation type %s not supported' % SIMULATION_TYPE)
        toc = time.time()-tic

        if model.solution is not None:

            with open(output_dir + '/' + output_file + '_flux.json', 'wb') as f:
                json.dump(model.get_metabolic_flux(), f)

            with open(output_dir + '/' + output_file + '_sol.pickle', 'wb') as f:
                cPickle.dump(model.solution, f)

            rows = []
            for rxn_id, v in model.solution.x_dict.items():
                # Save rxn, flux, keff (if rxn has keff)
                rxn = model.reactions.get_by_id(rxn_id)
                keff = rxn.keff if hasattr(rxn, 'keff') else np.nan
                rows.append({'rxn': rxn_id, 'v': v, 'keff': keff})
            df_result = pd.DataFrame(rows)
        else:
            df_result = pd.DataFrame([{'rxn': np.nan, 'v': np.nan,
                                       'keff': np.nan}])

        # ----------------------------------------------------
        # Write result to file
        #timestr = time.strftime("%Y%m%d")
        #filename = 'result_sampleme%d%s_%s_job%s_%d.csv' % (
        #    ME_PROTO_VER, simstr,
        #    timestr, SIMULATION_TYPE, sample)

        # Save solution
        #df_result.to_csv(filename, index=False)

        # Return
        result = {'basis': hs, 'xopt': xopt, 'muopt': muopt}
        return result, df_result, toc

    # Use correct indices and allow warm-start within a worker if needed
    sols = []

    hs = None
    for ind in inds:
        sol = simulate(ind, hs)
        res = sol[0]
        hs_new = res['basis']
        sols.append(sol)
        if hs_new is not None:
            hs = hs_new

    df_results = [sol[1] for sol in sols]
    times = [sol[2] for sol in sols]

    results = {'ind': inds, 'time': times}

    return results

# ============================================================
# Do work on data chunk
data = do_work(worker_tasks[rank])
print('Finished work by worker %d' % rank)

# ------------------------------------------------------------
# Gather results by root
data = comm.gather(data, root=0)

# ------------------------------------------------------------
# Report final result
if rank == 0:
    # Gather results
    print('Gathered data by root')

    # Save to pickle first
    timestr = time.strftime("%Y%m%d_%H%M")
    filename = 'results_sampleme%d%s_%s_job%s.pickle' % (
        ME_PROTO_VER, simstr,
        timestr, SIMULATION_TYPE)

    with open(filename, 'wb') as iofile:
        cPickle.dump(data, iofile)

    print('Saved gathered results to pickle file:', filename)

    # Save as dataframe
    inds = [i for res in data for i in res['ind']]
    times = [t for res in data for t in res['time']]

    df = pd.DataFrame({'ind': inds, 'time': times})

    filename_df = 'results_sampleme%d%s_job%s_%s.csv' % (
        ME_PROTO_VER, simstr,
        timestr, SIMULATION_TYPE)
    df.to_csv(filename_df, index=False)

    print('Saved results to csv:', filename_df)
