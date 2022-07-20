#!/usr/bin/env python3
from datetime import datetime
from os import chdir, getcwd, makedirs
from os.path import abspath, expanduser, isdir, isfile
from subprocess import check_output
from sys import argv, stdout
import argparse

# useful constants
VERSION = '0.0.1'
LOGFILE = None

# defaults
DEFAULT_PATH_ABM_HIV_COMMANDLINE = "/usr/local/bin/abm_hiv-HRSA_SD/abm_hiv_commandline.R"
DEFAULT_PATH_ABM_HIV_MODULES = "/usr/local/bin/abm_hiv-HRSA_SD/modules"
DEFAULT_FN_ABM_HIV_CALIBRATION = "abm_hiv_calibration_data.tsv"
DEFAULT_FN_ABM_HIV_DEMOGRAPHICS = "abm_hiv_demographic_data.tsv"
DEFAULT_FN_ABM_HIV_LOG = "log_abm_hiv.txt"
DEFAULT_FN_ABM_HIV_TIMES = "abm_hiv_times.tsv"
DEFAULT_FN_LOG = "log_favites.txt"

# check if user is just printing version
if '--version' in argv:
    print("FAVITES-COVID-Lite version %s" % VERSION); exit()

# return the current time as a string
def get_time():
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")

# print to the log (None implies stdout only)
def print_log(s='', end='\n'):
    tmp = "[%s] %s" % (get_time(), s)
    print(tmp, end=end); stdout.flush()
    if LOGFILE is not None:
        print(tmp, file=LOGFILE, end=end); LOGFILE.flush()

# run abm_hiv-HRSA_SD
def run_abm_hiv_hrsa_sd(outdir, abm_hiv_params_xlsx, abm_hiv_trans_start, abm_hiv_trans_end, abm_hiv_trans_time, path_abm_hiv_commandline=DEFAULT_PATH_ABM_HIV_COMMANDLINE, path_abm_hiv_modules=DEFAULT_PATH_ABM_HIV_MODULES, verbose=True):
    # run the R script
    command = ['Rscript', path_abm_hiv_commandline, path_abm_hiv_modules, abm_hiv_params_xlsx, str(abm_hiv_trans_start), str(abm_hiv_trans_end), str(abm_hiv_trans_time)]
    if verbose:
        print_log("Running abm_hiv-HRSA_SD Command: %s" % ' '.join(command))
    log_f = open('%s/%s' % (outdir,DEFAULT_FN_ABM_HIV_LOG), 'w')
    abm_out = check_output(command, stderr=log_f).decode(); log_f.flush()

    # separate the output
    if verbose:
        print_log("Parsing abm_hiv-HRSA_SD output...")
    log_text, abm_out = abm_out.split('[1] "Calibration metrics..."')
    log_f.write(log_text); log_f.close()
    calibration_data, abm_out = abm_out.strip().split('[1] "Transmission tree..."')
    transmission_data, abm_out = abm_out.strip().split('[1] "Sequence sample times..."')
    times_data, demographic_data = abm_out.strip().split('[1] "PLWH demographics..."')

    # parse the calibration data and write to file
    if verbose:
        print_log("Parsing calibration data...")
    calibration_data = [[v.strip() for v in l.strip().split()] for l in calibration_data.strip().splitlines()]
    fn = '%s/%s' % (outdir, DEFAULT_FN_ABM_HIV_CALIBRATION); f = open(fn, 'w')
    for row in calibration_data:
        f.write('\t'.join(row) + '\n')
    f.close()
    if verbose:
        print_log("Calibration data written to: %s" % fn)

    # parse the transmission data and write to file
    if verbose:
        print_log("Parsing transmission network...")
    fn = '%s/error_free_files/transmission_network.txt' % outdir; f = open(fn, 'w')
    for l in transmission_data.strip().splitlines()[1:]:
        f.write('\t'.join(v.strip() for v in l.strip().split()) + '\n')
    f.close()
    if verbose:
        print_log("Transmission network written to: %s" % fn)

    # parse the times data and write to file
    if verbose:
        print_log("Parsing all event times...")
    times_data = [[v.strip() for v in l.strip().split()] for l in times_data.strip().splitlines()]
    fn = '%s/%s' % (outdir, DEFAULT_FN_ABM_HIV_TIMES); f = open(fn, 'w')
    for row in times_data:
        f.write('\t'.join(row) + '\n')
    f.close()
    if verbose:
        print_log("All event times written to: %s" % fn)

    # parse the demographic data and write to file
    if verbose:
        print_log("Parsing demographic data...")
    fn = '%s/%s' % (outdir, DEFAULT_FN_ABM_HIV_DEMOGRAPHICS); f = open(fn, 'w')
    for l in demographic_data.strip().splitlines():
        f.write('\t'.join(v.strip() for v in l.strip().split()) + '\n')
    f.close()
    if verbose:
        print_log("Demographic data written to: %s" % fn)

# parse user args
def parse_args():
    # arg parser
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-o', '--output', required=True, type=str, help="Output Directory")
    parser.add_argument('--abm_hiv_params_xlsx', required=True, type=str, help="abm_hiv-HRSA-SD: Parameter XLSX File")
    parser.add_argument('--abm_hiv_trans_start', required=True, type=float, help="amb_hiv-HRSA-SD: Starting Transition Rate")
    parser.add_argument('--abm_hiv_trans_end', required=True, type=float, help="amb_hiv-HRSA-SD: Ending Transition Rate")
    parser.add_argument('--abm_hiv_trans_time', required=True, type=float, help="amb_hiv-HRSA-SD: Time (months) to Go from Starting to Ending Transition Rate")
    parser.add_argument('--sample_time_probs_csv', required=True, type=str, help="Sample Time Probabilities (CSV)")
    parser.add_argument('--gzip_output', action='store_true', help="Gzip Compress Output Files")
    parser.add_argument('--path_abm_hiv_commandline', required=False, type=str, default=DEFAULT_PATH_ABM_HIV_COMMANDLINE, help="Path to abm_hiv-HRSA_SD/abm_hiv_commandline.R")
    parser.add_argument('--path_abm_hiv_modules', required=False, type=str, default=DEFAULT_PATH_ABM_HIV_MODULES, help="Path to abm_hiv-HRSA_SD/modules")
    parser.add_argument('--version', action='store_true', help="Display Version")
    args = parser.parse_args()

    # fix/update args and return
    args.output = abspath(expanduser(args.output))
    args.abm_hiv_params_xlsx = abspath(expanduser(args.abm_hiv_params_xlsx))
    args.sample_time_probs_csv = abspath(expanduser(args.sample_time_probs_csv))
    args.path_abm_hiv_commandline = abspath(expanduser(args.path_abm_hiv_commandline))
    args.path_abm_hiv_modules = abspath(expanduser(args.path_abm_hiv_modules))
    return args

# check user args
def check_args(args):
    # check output directory 
    if isdir(args.output) or isfile(args.output):
        raise ValueError("Output directory exists: %s" % args.output)

    # check abm_hiv-HRSA_SD params XLSX
    if not isfile(args.abm_hiv_params_xlsx):
        raise ValueError("abm_hiv-HRSA_SD parameter XLSX file not found: %s" % args.abm_hiv_params_xlsx)

    # check starting transition rate
    if args.abm_hiv_trans_start <= 0:
        raise ValueError("abm_hiv-HRSA_SD starting transition rate must be positive: %s" % args.abm_hiv_trans_start)

    # check ending transition rate
    if args.abm_hiv_trans_end <= 0:
        raise ValueError("abm_hiv-HRSA_SD ending transition rate must be positive: %s" % args.abm_hiv_trans_end)

    # check time to go from starting to ending transition rate
    if args.abm_hiv_trans_time <= 0:
        raise ValueError("amb_hiv-HRSA_SD time from starting to ending transition rate must be positive: %s" % args.abm_hiv_trans_time)

    # check sample time probabilities
    if not isfile(args.sample_time_probs_csv):
        raise ValueError("Sample time probabilities (CSV) not found: %s" % args.sample_time_probs_csv)

    # check path_abm_hiv_commandline
    if not isfile(args.path_abm_hiv_commandline):
        raise ValueError("abm_hiv-HRSA_SD/abm_hiv_commandline.R not found: %s" % args.path_abm_hiv_commandline)

    # check path_abm_hiv_modules
    if not isdir(args.path_abm_hiv_modules):
        raise ValueError("abm_hiv-HRSA_SD/modules not found: %s" % args.path_abm_hiv_modules)

# load abm_hiv-HRSA_SD sample time probabilities
def load_sample_time_probs(sample_time_probs_fn):
    exit() # TODO

# main execution
if __name__ == "__main__":
    # parse/check user args and set things up
    args = parse_args(); check_args(args)
    makedirs(args.output, exist_ok=True)
    makedirs('%s/error_free_files' % args.output, exist_ok=True)
    LOGFILE_fn = "%s/%s" % (args.output, DEFAULT_FN_LOG)
    LOGFILE = open(LOGFILE_fn, 'w')

    # print run info to log
    print_log("===== RUN INFORMATION =====")
    print_log("=== Command Information ===")
    print_log("FAVITES-ABM-HIV-SD-Lite Version: %s" % VERSION)
    print_log("FAVITES-ABM-HIV-SD-Lite Command: %s" % ' '.join(argv))
    print_log()
    print_log("=== Command Arguments ===")
    print_log("Output Directory: %s" % args.output)
    print_log("Gzip Output: %s" % args.gzip_output)
    print_log(); print_log()

    # simulate contact and transmission network
    print_log("===== CONTACT AND TRANSMISSION NETWORK =====")
    print_log("=== Contact and Transmission Network Arguments ===")
    print_log("abm_hiv-HRSA_SD: Parameter XLSX: %s" % args.abm_hiv_params_xlsx)
    print_log("abm_hiv-HRSA_SD: Starting Transition Rate: %s" % args.abm_hiv_trans_start)
    print_log("abm_hiv-HRSA_SD: Ending Transition Rate: %s" % args.abm_hiv_trans_end)
    print_log("abm_hiv-HRSA_SD: Time (months) to Go from Starting to Ending Transition Rate: %s" % args.abm_hiv_trans_time)
    print_log("Sample Time Probabilities (CSV): %s" % args.sample_time_probs_csv)
    print_log()
    print_log("=== abm_hiv-HRSA_SD Progress ===")
    abm_out = run_abm_hiv_hrsa_sd(
        args.output,                   # FAVITES-ABM-HIV-SD-Lite output directory
        args.abm_hiv_params_xlsx,      # parameter XLSX file
        args.abm_hiv_trans_start,      # starting transition rate
        args.abm_hiv_trans_end,        # ending transition rate
        args.abm_hiv_trans_time,       # time to go from starting to ending transition rate
        args.path_abm_hiv_commandline, # path to abm_hiv-HRSA_SD/abm_hiv_commandline.R script
        args.path_abm_hiv_modules)     # path to abm_hiv-HRSA_SD/modules folder
    print_log()
