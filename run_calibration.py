#!/usr/bin/env python3
from run_favites_lite import DEFAULT_PATH_ABM_HIV_COMMANDLINE, DEFAULT_PATH_ABM_HIV_MODULES, DEFAULT_PATH_COATRAN_CONSTANT
from csv import reader
from datetime import datetime
from glob import glob
from openpyxl import load_workbook
from os import chdir, cpu_count, getcwd, makedirs
from os.path import abspath, expanduser, isdir, isfile
from random import randint, seed
from scipy.optimize import minimize
from shutil import make_archive, rmtree
from subprocess import check_output
from sys import argv, stdout
from zipfile import ZIP_DEFLATED, ZipFile
import argparse

# useful constants
LOGFILE = None
DEFAULT_FN_LOG = "log_favites_calibration.txt"
MAX_RNG_SEED = 2147483647
SCIPY_MINIMIZE_METHODS = {'Nelder-Mead', 'L-BFGS-B', 'TNC', 'SLSQP', 'Powell', 'trust-constr', 'COBYLA', 'COBYQA'}

# return the current time as a string
def get_time():
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")

# print to the log (None implies stdout only)
def print_log(s='', end='\n'):
    tmp = "[%s] %s" % (get_time(), s)
    print(tmp, end=end); stdout.flush()
    if LOGFILE is not None:
        print(tmp, file=LOGFILE, end=end); LOGFILE.flush()

# parse user args
def parse_args():
    # arg parser
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-o', '--output', required=True, type=str, help="Output Directory")
    parser.add_argument('--calibration_mode', required=True, type=str, help="Calibration Mode (epi, epi+genetic)")
    parser.add_argument('--calibration_csv', required=True, type=str, help="Calibration CSV")
    parser.add_argument('--sim_start_time', required=True, type=float, help="Simulation Start Time (year)")
    parser.add_argument('--abm_hiv_params_xlsx', required=True, type=str, help="abm_hiv-HRSA_SD: Parameter XLSX File")
    parser.add_argument('--abm_hiv_sd_demographics_csv', required=True, type=str, help="abm_hiv-HRSA_SD: Demographics CSV File")
    parser.add_argument('--abm_hiv_trans_start', required=True, type=float, help="abm_hiv-HRSA_SD: Starting Transition Rate")
    parser.add_argument('--abm_hiv_trans_end', required=True, type=float, help="abm_hiv-HRSA_SD: Ending Transition Rate")
    parser.add_argument('--abm_hiv_trans_time', required=True, type=float, help="abm_hiv-HRSA_SD: Time (months) to Go from Starting to Ending Transition Rate")
    parser.add_argument('--sample_time_probs_csv', required=True, type=str, help="Sample Time Probabilities (CSV)")
    parser.add_argument('--coatran_eff_pop_size', required=True, type=float, help="CoaTran (Constant): Effective Population Size")
    parser.add_argument('--mutation_tree_seed', required=True, type=str, help="Mutation Tree: Seed (Newick/Nexus File)")
    parser.add_argument('--time_tree_seed', required=True, type=str, help="Time Tree: Seed (Newick/Nexus File)")
    parser.add_argument('--time_tree_tmrca', required=True, type=float, help="Time Tree: Time of Most Recent Common Ancestor (tMRCA; time of root of seed time tree; year)")
    parser.add_argument('--time_tree_only_include_mapped', action='store_true', help="Time Tree: Only Include Individuals in ABM Map in Seed Tree")
    parser.add_argument('--scipy_minimize_method', required=False, type=str, default='Powell', help="SciPy Minimize Optimization Method (options: %s)" % ', '.join(sorted(SCIPY_MINIMIZE_METHODS)))
    parser.add_argument('--num_reps_per_score', required=False, type=int, default=5, help="Number of replicates to run per scoring of a parameter set")
    parser.add_argument('--max_num_threads', required=False, type=int, default=cpu_count(), help="Max number of threads to use")
    parser.add_argument('--zip_output', action='store_true', help="Gzip Compress Output Files")
    parser.add_argument('--path_abm_hiv_commandline', required=False, type=str, default=DEFAULT_PATH_ABM_HIV_COMMANDLINE, help="Path to abm_hiv-HRSA_SD/abm_hiv_commandline.R")
    parser.add_argument('--path_abm_hiv_modules', required=False, type=str, default=DEFAULT_PATH_ABM_HIV_MODULES, help="Path to abm_hiv-HRSA_SD/modules")
    parser.add_argument('--path_coatran_constant', required=False, type=str, default=DEFAULT_PATH_COATRAN_CONSTANT, help="Path to coatran_constant")
    parser.add_argument('--version', action='store_true', help="Display Version")
    args = parser.parse_args()

    # fix/update args
    args.output = abspath(expanduser(args.output))
    args.calibration_csv = abspath(expanduser(args.calibration_csv))
    args.abm_hiv_params_xlsx = abspath(expanduser(args.abm_hiv_params_xlsx))
    args.abm_hiv_sd_demographics_csv = abspath(expanduser(args.abm_hiv_sd_demographics_csv))
    args.sample_time_probs_csv = abspath(expanduser(args.sample_time_probs_csv))
    args.time_tree_seed = abspath(expanduser(args.time_tree_seed))
    args.path_abm_hiv_commandline = abspath(expanduser(args.path_abm_hiv_commandline))
    args.path_abm_hiv_modules = abspath(expanduser(args.path_abm_hiv_modules))
    args.scipy_minimize_method = args.scipy_minimize_method.strip()

    # check input files and return
    for fn in [args.calibration_csv, args.abm_hiv_params_xlsx, args.abm_hiv_sd_demographics_csv, args.sample_time_probs_csv, args.time_tree_seed]:
        if not isfile(fn):
            raise ValueError("File not found: %s" % fn)
    return args

# check user args (most will be checked by run_favites_lite.py)
def check_args(args):
    # check GNU parallel
    if not check_output(['parallel', '-h']).decode().strip().startswith('Usage:'):
        raise RuntimeError("GNU `parallel` does not seem to be installed.\nInstall with e.g.: sudo apt-get install -y parallel")

    # check output directory 
    if isdir(args.output) or isfile(args.output):
        raise ValueError("Output directory exists: %s" % args.output)

    # check calibration CSV
    if not isfile(args.calibration_csv):
        raise ValueError("Calibration CSV not found: %s" % args.calibration_csv)

    # check scipy.minimize optimization method
    if args.scipy_minimize_method not in SCIPY_MINIMIZE_METHODS:
        raise ValueError("Invalid scipy.minimize method (%s). Options: %s" % (args.scipy_minimize_method, ', '.join(sorted(SCIPY_MINIMIZE_METHODS))))

    # check num reps per score
    if args.num_reps_per_score < 1:
        raise ValueError("Number of replicates per scoring must be at least 1: %s" % args.num_reps_per_score)

    # check max num threads
    if args.max_num_threads < 1 or args.max_num_threads > cpu_count():
        raise ValueError("Max number of threads must be at least 1 and at most the total number of CPUs (%d): %s" % (cpu_count(), args.max_num_threads))

# run calibration
def run_calibration(
        calibration_csv, output, sim_start_time,
        abm_hiv_params_xlsx, abm_hiv_sd_demographics_csv, abm_hiv_trans_start, abm_hiv_trans_end, abm_hiv_trans_time,
        sample_time_probs_csv, coatran_eff_pop_size, time_tree_seed, time_tree_tmrca, time_tree_only_include_mapped,
        mutation_tree_seed,
        scipy_minimize_method, # options in `SCIPY_MINIMIZE_METHODS`
        num_reps_per_score, max_num_threads,
        calibration_mode, # "epi" or "epi+genetic"
        path_abm_hiv_commandline=DEFAULT_PATH_ABM_HIV_COMMANDLINE, path_abm_hiv_modules=DEFAULT_PATH_ABM_HIV_MODULES, path_coatran_constant=DEFAULT_PATH_COATRAN_CONSTANT,
        zip_output=False,
    ):
    # prep calibration parameters
    calibration_mode_parts = set(calibration_mode.lower().split('+'))
    x0 = list(); bounds = list(); override_cli_args = list()
    if 'epi' in calibration_mode_parts:
        # IDU injections per month
        x0.append(0.05)
        bounds.append((0, 1))
        override_cli_args.append('--abm_xlsx_idu_injections_per_month')

        # heterosexual partners per person per month
        x0.append(1.1)
        bounds.append((0.8, 1.5))
        override_cli_args.append('--abm_xlsx_het_partners_per_person_per_month')

        # heterosexual assortativity (HIV status)
        x0.append(0.7)
        bounds.append((0, 1))
        override_cli_args.append('--abm_xlsx_het_assortativity_hiv_status')

        # MSM partners per person per month
        x0.append(1.1)
        bounds.append((0.8, 1.5))
        override_cli_args.append('--abm_xlsx_msm_partners_per_person_per_month')

        # MSM assortativity (HIV status)
        x0.append(0.7)
        bounds.append((0, 1))
        override_cli_args.append('--abm_xlsx_msm_assortativity_hiv_status')

        # MSM HIV monthly transmission probability between partners multiplier
        x0.append(1)
        bounds.append((1, 2))
        override_cli_args.append('--abm_xlsx_msm_hiv_monthly_transmission_prob_mult')

        if 'genetic' in calibration_mode_parts:
            # heterosexual assortativity (risk)
            x0.append(0.2)
            bounds.append((0, 1))
            override_cli_args.append('--abm_xlsx_het_assortativity_risk')

            # MSM assortativity (risk)
            x0.append(0.2)
            bounds.append((0, 1))
            override_cli_args.append('--abm_xlsx_msm_assortativity_risk')

            # MSMW percentage of MSM who are MSMW
            x0.append(0)
            bounds.append((0, 0.2))
            override_cli_args.append('--abm_xlsx_msmw_percentage_msm_msmw')

    # check for validity
    if len(x0) == 0:
        raise ValueError("Invalid calibration mode: %s" % calibration_mode)

    # prep calibration execution
    iter_num = 0
    run_favites_lite_path = '%s/%s' % ('/'.join(__file__.replace('/./','/').split('/')[:-1]), 'run_favites_lite.py')
    score_simulation_output_path = '%s/%s' % ('/'.join(__file__.replace('/./','/').split('/')[:-1]), 'score_simulation_output.py')

    # nested optimization function
    def opt_func(x):
        # prep FAVITES run
        nonlocal iter_num; iter_num += 1
        curr_outdir = '%s/favites.output.%s' % (output, str(iter_num).zfill(3))
        makedirs(curr_outdir, exist_ok=True)
        rep_nums = [str(i).zfill(3) for i in range(1, num_reps_per_score+1)]
        seed(); rng_seed_base = randint(0, 99999)
        run_favites_lite_command = [
            'parallel', '--jobs', str(max_num_threads), '--link',
            run_favites_lite_path,
            '--output', '%s/{1}' % curr_outdir,
            '--sim_start_time', str(sim_start_time),
            '--abm_hiv_params_xlsx', abm_hiv_params_xlsx,
            '--abm_hiv_sd_demographics_csv', abm_hiv_sd_demographics_csv,
            '--abm_hiv_trans_start', str(abm_hiv_trans_start),
            '--abm_hiv_trans_end', str(abm_hiv_trans_end),
            '--abm_hiv_trans_time', str(abm_hiv_trans_time),
            '--sample_time_probs_csv', sample_time_probs_csv,
            '--coatran_eff_pop_size', str(coatran_eff_pop_size),
            '--mutation_tree_seed', mutation_tree_seed,
            '--time_tree_seed', time_tree_seed,
            '--time_tree_tmrca', str(time_tree_tmrca),
            '--path_abm_hiv_commandline', path_abm_hiv_commandline,
            '--path_abm_hiv_modules', path_abm_hiv_modules,
            '--path_coatran_constant', path_coatran_constant,
            '--abm_xlsx_rng_seed', '{2}',
        ]
        if time_tree_only_include_mapped:
            run_favites_lite_command.append('--time_tree_only_include_mapped')

        # add calibration parameters
        for i in range(len(override_cli_args)):
            run_favites_lite_command += [override_cli_args[i], str(x[i])]

        # run FAVITES
        run_favites_lite_command += [':::'] + rep_nums
        run_favites_lite_command += [':::'] + [str(rng_seed_base+int(v)) for v in rep_nums]
        print_log("Running FAVITES iteration %d: %s" % (iter_num, ' '.join(run_favites_lite_command)))
        try:
            o = check_output(run_favites_lite_command)
        except:
            print_log("WARNING: At least one run crashed in FAVITES iteration %d" % iter_num)

        # calculate optimization function score
        score_simulation_output_command = [
            'parallel', '--jobs', str(max_num_threads),
            score_simulation_output_path,
            '-s', '%s/{1}' % curr_outdir,
            '-c', calibration_csv,
            '-m', calibration_mode,
            '-o', '%s/{1}/score.txt' % curr_outdir,
        ] + [':::'] + rep_nums
        print_log("Calculating optimization function score: %s" % ' '.join(score_simulation_output_command))
        try:
            check_output(score_simulation_output_command)
        except:
            print_log("WARNING: Failed to score at least one run in FAVITES iteration %d" % iter_num)
        if len(list(glob('%s/*/score.txt' % curr_outdir))) == 0:
            error_message = "All simulation replicates crashed in FAVITES iteration %d" % iter_num
            print_log(error_message)
            raise RuntimeError(error_message)
        score = min(float(open(fn).read().split()[-1]) for fn in glob('%s/*/score.txt' % curr_outdir))
        print_log("FAVITES iteration %d average score: %s" % (iter_num, score))
        if zip_output:
            print_log("Zipping output..."); make_archive(curr_outdir, 'zip', curr_outdir); rmtree(curr_outdir)
        return score

    # run calibration
    minimize(opt_func, x0, bounds=bounds, method=scipy_minimize_method)

# main execution
if __name__ == "__main__":
    # parse/check user args and set things up
    args = parse_args(); check_args(args)
    makedirs(args.output, exist_ok=True)
    LOGFILE_fn = "%s/%s" % (args.output, DEFAULT_FN_LOG)
    LOGFILE = open(LOGFILE_fn, 'w')

    # print run info to log
    print_log("===== RUN INFORMATION =====")
    print_log("Calibration command: %s" % ' '.join(argv))
    print_log(); print_log()

    # run calibration
    print_log("===== CALIBRATION =====")
    print_log("Running calibration...")
    run_calibration(
        args.calibration_csv, args.output, args.sim_start_time,
        args.abm_hiv_params_xlsx, args.abm_hiv_sd_demographics_csv, args.abm_hiv_trans_start, args.abm_hiv_trans_end, args.abm_hiv_trans_time,
        args.sample_time_probs_csv, args.coatran_eff_pop_size, args.time_tree_seed, args.time_tree_tmrca, args.time_tree_only_include_mapped,
        args.mutation_tree_seed,
        args.scipy_minimize_method,
        args.num_reps_per_score, args.max_num_threads,
        args.calibration_mode.lower().strip(),
        args.path_abm_hiv_commandline, args.path_abm_hiv_modules, args.path_coatran_constant,
        zip_output=args.zip_output,
    )
