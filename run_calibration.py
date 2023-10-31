#!/usr/bin/env python3

# TODO REWRITE TO CALL EXTERNAL "SIMULATION FOLDER SCORING" HELPER SCRIPT INTEAD OF SCORING INTERNALLY

from run_favites_lite import DEFAULT_PATH_ABM_HIV_COMMANDLINE, DEFAULT_PATH_ABM_HIV_MODULES, DEFAULT_PATH_COATRAN_CONSTANT
from csv import reader
from datetime import datetime
from openpyxl import load_workbook
from os import chdir, getcwd, makedirs
from os.path import abspath, expanduser, isdir, isfile
from random import randint
from scipy.optimize import minimize
from subprocess import check_output
from sys import argv, stdout
import argparse

# useful constants
LOGFILE = None
DEFAULT_FN_LOG = "log_favites_calibration.txt"
MAX_RNG_SEED = 2147483647

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
    parser.add_argument('--calibration_csv', required=True, type=str, help="Calibration CSV")
    parser.add_argument('--sim_start_time', required=True, type=float, help="Simulation Start Time (year)")
    parser.add_argument('--abm_hiv_params_xlsx', required=True, type=str, help="abm_hiv-HRSA_SD: Parameter XLSX File")
    parser.add_argument('--abm_hiv_sd_demographics_csv', required=True, type=str, help="abm_hiv-HRSA_SD: Demographics CSV File")
    parser.add_argument('--abm_hiv_trans_start', required=True, type=float, help="abm_hiv-HRSA_SD: Starting Transition Rate")
    parser.add_argument('--abm_hiv_trans_end', required=True, type=float, help="abm_hiv-HRSA_SD: Ending Transition Rate")
    parser.add_argument('--abm_hiv_trans_time', required=True, type=float, help="abm_hiv-HRSA_SD: Time (months) to Go from Starting to Ending Transition Rate")
    parser.add_argument('--sample_time_probs_csv', required=True, type=str, help="Sample Time Probabilities (CSV)")
    parser.add_argument('--coatran_eff_pop_size', required=True, type=float, help="CoaTran (Constant): Effective Population Size")
    parser.add_argument('--time_tree_seed', required=True, type=str, help="Time Tree: Seed (Newick File)")
    parser.add_argument('--time_tree_tmrca', required=True, type=float, help="Time Tree: Time of Most Recent Common Ancestor (tMRCA; time of root of seed time tree; year)")
    parser.add_argument('--time_tree_only_include_mapped', action='store_true', help="Time Tree: Only Include Individuals in ABM Map in Seed Tree")
    parser.add_argument('--mutation_rate_loc', required=True, type=float, help="Mutation Rate: Truncated Normal Location (mutations/month)")
    parser.add_argument('--mutation_rate_scale', required=True, type=float, help="Mutation Rate: Truncated Normal Scale (mutations/month)")
    parser.add_argument('--gzip_output', action='store_true', help="Gzip Compress Output Files")
    parser.add_argument('--path_abm_hiv_commandline', required=False, type=str, default=DEFAULT_PATH_ABM_HIV_COMMANDLINE, help="Path to abm_hiv-HRSA_SD/abm_hiv_commandline.R")
    parser.add_argument('--path_abm_hiv_modules', required=False, type=str, default=DEFAULT_PATH_ABM_HIV_MODULES, help="Path to abm_hiv-HRSA_SD/modules")
    parser.add_argument('--path_coatran_constant', required=False, type=str, default=DEFAULT_PATH_COATRAN_CONSTANT, help="Path to coatran_constant")
    parser.add_argument('--version', action='store_true', help="Display Version")
    args = parser.parse_args()

    # fix/update args and return
    args.output = abspath(expanduser(args.output))
    args.calibration_csv = abspath(expanduser(args.calibration_csv))
    args.abm_hiv_params_xlsx = abspath(expanduser(args.abm_hiv_params_xlsx))
    args.abm_hiv_sd_demographics_csv = abspath(expanduser(args.abm_hiv_sd_demographics_csv))
    args.sample_time_probs_csv = abspath(expanduser(args.sample_time_probs_csv))
    args.time_tree_seed = abspath(expanduser(args.time_tree_seed))
    args.path_abm_hiv_commandline = abspath(expanduser(args.path_abm_hiv_commandline))
    args.path_abm_hiv_modules = abspath(expanduser(args.path_abm_hiv_modules))
    return args

# check user args (most will be checked by run_favites_lite.py)
def check_args(args):
    # check output directory 
    if isdir(args.output) or isfile(args.output):
        raise ValueError("Output directory exists: %s" % args.output)

    # check calibration CSV
    if not isfile(args.calibration_csv):
        raise ValueError("Calibration CSV not found: %s" % args.calibration_csv)

# load calibration data from CSV; calibration_data[description] = (value, weight)
def load_calibration_data(calibration_csv):
    calibration_data = dict()
    for row in reader(open(calibration_csv)):
        d, v, w = [x.strip() for x in row]
        if d.lower() == 'description':
            continue # skip header row
        if d in calibration_data:
            raise ValueError("Duplicate entry in calibration CSV: %s" % d)
        calibration_data[d] = (float(v), float(w))
    return calibration_data

# load calibration data output by the ABM R code; data[metric][year][subgroup] = count
def load_abm_calibration_output(abm_calibration_tsv):
    data = dict()
    for row in reader(open(abm_calibration_tsv), delimiter='\t'):
        metric, month, subgroup, stat = [x.strip() for x in row]
        if metric.lower() == 'metric':
            continue # skip header row
        if metric not in data:
            data[metric] = dict()
        year = int(month) // 12
        if year not in data[metric]:
            data[metric][year] = dict()
        if subgroup not in data[metric][year]:
            data[metric][year][subgroup] = 0
        data[metric][year][subgroup] += int(stat)
    return data

# run calibration
def run_calibration(
        calibration_data, output, sim_start_time,
        abm_hiv_params_xlsx, abm_hiv_sd_demographics_csv, abm_hiv_trans_start, abm_hiv_trans_end, abm_hiv_trans_time,
        sample_time_probs_csv, coatran_eff_pop_size, time_tree_seed, time_tree_tmrca, time_tree_only_include_mapped,
        mutation_rate_loc, mutation_rate_scale,
        path_abm_hiv_commandline=DEFAULT_PATH_ABM_HIV_COMMANDLINE, path_abm_hiv_modules=DEFAULT_PATH_ABM_HIV_MODULES, path_coatran_constant=DEFAULT_PATH_COATRAN_CONSTANT,
    ):
    # prep calibration
    x0 = [0] # TODO REPLACE WITH ACTUAL STARTING VALUES BASED ON calibration_data
    bounds = [
        (-100, 100), # TODO SET BOUNDS FOR x0[0], x0[1], etc.
    ]
    options = {
        'maxiter': 2, # max number of iterations (TODO change from 1 to some other value; maybe function parameter with default value?)
    }
    iter_num = 0
    run_favites_lite_path = '%s/%s' % ('/'.join(__file__.replace('/./','/').split('/')[:-1]), 'run_favites_lite.py')
    abm_xlsx_dir = '%s/abm_xlsx' % output
    makedirs(abm_xlsx_dir, exist_ok=True)

    # nested optimization function
    def opt_func(x):
        # set up FAVITES run
        nonlocal iter_num; iter_num += 1
        curr_outdir = '%s/favites.output.%s' % (output, str(iter_num).zfill(3))
        curr_data_xlsx = '%s/abm.data.%s.xlsx' % (abm_xlsx_dir, str(iter_num).zfill(3))
        wb = load_workbook(abm_hiv_params_xlsx) # load original ABM XLSX
        wb['High Level Pop + Sim Features']['B4'] = randint(0, MAX_RNG_SEED) # replace ABM seed
        wb['High Level Pop + Sim Features']['B4'].number_format = '0' # ABM seed is int, not string
        wb.save(curr_data_xlsx)

        # run FAVITES
        run_favites_lite_command = [
            run_favites_lite_path,
            '--output', curr_outdir,
            '--sim_start_time', str(sim_start_time),
            '--abm_hiv_params_xlsx', curr_data_xlsx,
            '--abm_hiv_sd_demographics_csv', abm_hiv_sd_demographics_csv,
            '--abm_hiv_trans_start', str(abm_hiv_trans_start),
            '--abm_hiv_trans_end', str(abm_hiv_trans_end),
            '--abm_hiv_trans_time', str(abm_hiv_trans_time),
            '--sample_time_probs_csv', sample_time_probs_csv,
            '--coatran_eff_pop_size', str(coatran_eff_pop_size),
            '--time_tree_seed', time_tree_seed,
            '--time_tree_tmrca', str(time_tree_tmrca),
            '--mutation_rate_loc', str(mutation_rate_loc),
            '--mutation_rate_scale', str(mutation_rate_scale),
            '--path_abm_hiv_commandline', path_abm_hiv_commandline,
            '--path_abm_hiv_modules', path_abm_hiv_modules,
            '--path_coatran_constant', path_coatran_constant,
        ]
        if time_tree_only_include_mapped:
            run_favites_lite_command.append('--time_tree_only_include_mapped')
        print_log("Running FAVITES iteration %d: %s" % (iter_num, ' '.join(run_favites_lite_command)))
        o = check_output(run_favites_lite_command)

        # calculate optimization function score TODO REFACTOR INTO OWN FUNCTION
        abm_calibration_data = load_abm_calibration_output('%s/abm_hiv_calibration_data.tsv' % curr_outdir)
        score = 0
        for cal_key, cal_tup in calibration_data.items():
            sim_val = None; cal_val, cal_w = cal_tup
            if '_' in cal_key:
                try:
                    year, risk = cal_key.split('_')
                    year = int(year)
                    risk = risk.replace(' & ', 'and')
                    if risk.lower() == 'other':
                        risk = 'other' # 'other' is lowercase in the ABM R code output
                    sim_val = abm_calibration_data['newinfects_agg'][year-sim_start_time][risk]
                except:
                    pass # sim_val will remain None, so ValueError below will be thrown
            if sim_val is None:
                raise ValueError("Unknown calibration key: %s" % cal_key)
            else:
                score += cal_w * ((sim_val-cal_val)**2)
        score = score**0.5
        print_log("FAVITES iteration %d score: %s" % (iter_num, score))
        return score

    # run calibration
    minimize(opt_func, x0, bounds=bounds, options=options)

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
    print_log("Loading calibration data: %s" % args.calibration_csv)
    calibration_data = load_calibration_data(args.calibration_csv)
    print_log("Running calibration...")
    run_calibration(
        calibration_data, args.output, args.sim_start_time,
        args.abm_hiv_params_xlsx, args.abm_hiv_sd_demographics_csv, args.abm_hiv_trans_start, args.abm_hiv_trans_end, args.abm_hiv_trans_time,
        args.sample_time_probs_csv, args.coatran_eff_pop_size, args.time_tree_seed, args.time_tree_tmrca, args.time_tree_only_include_mapped,
        args.mutation_rate_loc, args.mutation_rate_scale,
        args.path_abm_hiv_commandline, args.path_abm_hiv_modules, args.path_coatran_constant,
    )
    pass # TODO CONTINUE HERE
