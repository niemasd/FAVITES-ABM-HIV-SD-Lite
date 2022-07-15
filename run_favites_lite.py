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
def run_abm_hiv_hrsa_sd(outdir, abm_hiv_params_xlsx, abm_hiv_trans_start, abm_hiv_trans_end, abm_hiv_trans_time, path_abm_hiv_commandline=DEFAULT_PATH_ABM_HIV_COMMANDLINE, path_abm_hiv_modules=DEFAULT_PATH_ABM_HIV_MODULES):
    # run the R script
    command = ['Rscript', path_abm_hiv_commandline, path_abm_hiv_modules, abm_hiv_params_xlsx, str(abm_hiv_trans_start), str(abm_hiv_trans_end), str(abm_hiv_trans_time)]
    print_log("abm_hiv-HRSA_SD Command: %s" % ' '.join(command))
    log_f = open('%s/%s' % (outdir,DEFAULT_FN_ABM_HIV_LOG), 'w')
    abm_out = check_output(command, stderr=log_f).decode(); log_f.flush()

    # separate the output
    log_text, abm_out = abm_out.split('[1] "Calibration metrics..."')
    log_f.write(log_text); log_f.close()
    calibration_data, abm_out = abm_out.strip().split('[1] "Transmission tree..."')
    transmission_data, times_data = abm_out.strip().split('[1] "Sequence sample times..."')

    # parse the calibration data and write to file
    calibration_data = [[v.strip() for v in l.strip().split()] for l in calibration_data.strip().splitlines()[1:]]
    f = open('%s/%s' % (outdir, DEFAULT_FN_ABM_HIV_CALIBRATION), 'w')
    f.write('\t'.join(calibration_data[0]) + '\n')
    for row in calibration_data[2:]:
        f.write('\t'.join(row[1:]) + '\n')
    f.close()

    # parse the transmission data and write to file
    f = open('%s/error_free_files/transmission_network.txt' % outdir, 'w')
    for l in transmission_data.strip().splitlines()[3:]:
        f.write('\t'.join(v.strip() for v in l.strip().split()[1:]) + '\n')
    f.close()

    # parse the times data and write to file
    times_data = [[v.strip() for v in l.strip().split()] for l in times_data.strip().splitlines()[1:]]
    f = open('%s/%s' % (outdir, DEFAULT_FN_ABM_HIV_TIMES), 'w')
    f.write('\t'.join(times_data[0]) + '\n')
    for row in times_data[2:]:
        f.write('\t'.join(row[1:]) + '\n')
    f.close()

# parse user args
def parse_args():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-o', '--output', required=True, type=str, help="Output Directory")
    parser.add_argument('--abm_hiv_params_xlsx', required=True, type=str, help="abm_hiv-HRSA-SD: Parameter XLSX File")
    parser.add_argument('--abm_hiv_trans_start', required=True, type=float, help="amb_hiv-HRSA-SD: Starting Transition Rate")
    parser.add_argument('--abm_hiv_trans_end', required=True, type=float, help="amb_hiv-HRSA-SD: Ending Transition Rate")
    parser.add_argument('--abm_hiv_trans_time', required=True, type=float, help="amb_hiv-HRSA-SD: Time (months) to Go from Starting to Ending Transition Rate")
    parser.add_argument('--gzip_output', action='store_true', help="Gzip Compress Output Files")
    parser.add_argument('--path_abm_hiv_commandline', required=False, type=str, default=DEFAULT_PATH_ABM_HIV_COMMANDLINE, help="Path to abm_hiv-HRSA_SD/abm_hiv_commandline.R")
    parser.add_argument('--path_abm_hiv_modules', required=False, type=str, default=DEFAULT_PATH_ABM_HIV_MODULES, help="Path to abm_hiv-HRSA_SD/modules")
    parser.add_argument('--version', action='store_true', help="Display Version")
    return parser.parse_args()

# main execution
if __name__ == "__main__":
    # parse and check user args
    args = parse_args()
    if isdir(args.output) or isfile(args.output):
        raise ValueError("Output directory exists: %s" % args.output)
    args.output = abspath(expanduser(args.output))

    # set up output directory
    makedirs(args.output, exist_ok=True)
    makedirs('%s/error_free_files' % args.output, exist_ok=True)
    LOGFILE_fn = "%s/%s" % (args.output, DEFAULT_FN_LOG)
    LOGFILE = open(LOGFILE_fn, 'w')

    # print run info to log
    print_log("===== RUN INFORMATION =====")
    print_log("FAVITES-ABM-HIV-SD-Lite Version: %s" % VERSION)
    print_log("FAVITES-ABM-HIV-SD-Lite Command: %s" % ' '.join(argv))
    print_log("Output Directory: %s" % args.output)
    print_log()

    # simulate contact and transmission network
    print_log("===== CONTACT AND TRANSMISSION NETWORK =====")
    abm_out = run_abm_hiv_hrsa_sd(
        args.output,                   # FAVITES-ABM-HIV-SD-Lite output directory
        args.abm_hiv_params_xlsx,      # parameter XLSX file
        args.abm_hiv_trans_start,      # starting transition rate
        args.abm_hiv_trans_end,        # ending transition rate
        args.abm_hiv_trans_time,       # time to go from starting to ending transition rate
        args.path_abm_hiv_commandline, # path to abm_hiv-HRSA_SD/abm_hiv_commandline.R script
        args.path_abm_hiv_modules)     # path to abm_hiv-HRSA_SD/modules folder
    print_log()
