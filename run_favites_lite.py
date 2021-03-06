#!/usr/bin/env python3
from datetime import datetime
from math import exp
from os import chdir, getcwd, makedirs
from os.path import abspath, expanduser, isdir, isfile
from random import random
from scipy.stats import truncnorm
from subprocess import check_output
from sys import argv, stdout
from treesap import nonhomogeneous_yule_tree
from treeswift import read_tree_newick
import argparse

# useful constants
VERSION = '0.0.1'
LOGFILE = None
SAMPLE_TIME_PROB_COLUMNS = ['gender', 'risk', 'race', 'agerange', 'timetype', 'probability']
DEMOGRAPHICS_COLUMNS = ['id', 'gender', 'risk', 'age', 'race']
AGE_RANGES = ['13-24', '25-54', '55-100']
AGE_RANGES_FLOAT = [(float(v.split('-')[0]), float(v.split('-')[1])) for v in AGE_RANGES]

# defaults
DEFAULT_PATH_ABM_HIV_COMMANDLINE = "/usr/local/bin/abm_hiv-HRSA_SD/abm_hiv_commandline.R"
DEFAULT_PATH_ABM_HIV_MODULES = "/usr/local/bin/abm_hiv-HRSA_SD/modules"
DEFAULT_PATH_COATRAN_CONSTANT = "coatran_constant"
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
    end_time = float(log_text.split('"Month:')[-1].split('"')[0])
    if verbose:
        print_log("Simulation end time (months): %s" % end_time)
    
    # parse the calibration data and write to file
    if verbose:
        print_log("Parsing calibration data...")
    calibration_data = [[v.strip() for v in l.strip().split()] for l in calibration_data.strip().splitlines()]
    calibration_fn = '%s/%s' % (outdir, DEFAULT_FN_ABM_HIV_CALIBRATION); f = open(calibration_fn, 'w')
    for row in calibration_data:
        f.write('\t'.join(row) + '\n')
    f.close()
    if verbose:
        print_log("Calibration data written to: %s" % calibration_fn)

    # parse the transmission data and write to file
    if verbose:
        print_log("Parsing transmission network...")
    transmission_fn = '%s/error_free_files/transmission_network.txt' % outdir; f = open(transmission_fn, 'w')
    for l in transmission_data.strip().splitlines()[1:]:
        f.write('\t'.join(v.strip() for v in l.strip().split()) + '\n')
    f.close()
    if verbose:
        print_log("Transmission network written to: %s" % transmission_fn)

    # parse the times data and write to file
    if verbose:
        print_log("Parsing all event times...")
    times_data = [[v.strip() for v in l.strip().split()] for l in times_data.strip().splitlines()]
    all_times_fn = '%s/%s' % (outdir, DEFAULT_FN_ABM_HIV_TIMES); f = open(all_times_fn, 'w')
    for row in times_data:
        f.write('\t'.join(row) + '\n')
    f.close()
    if verbose:
        print_log("All event times written to: %s" % all_times_fn)

    # parse the demographic data and write to file
    if verbose:
        print_log("Parsing demographic data...")
    demographic_fn = '%s/%s' % (outdir, DEFAULT_FN_ABM_HIV_DEMOGRAPHICS); f = open(demographic_fn, 'w')
    for l in demographic_data.strip().splitlines():
        f.write('\t'.join(v.strip() for v in l.strip().split()) + '\n')
    f.close()
    if verbose:
        print_log("Demographic data written to: %s" % demographic_fn)

    # return output filenames
    return end_time, calibration_fn, transmission_fn, all_times_fn, demographic_fn

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
    parser.add_argument('--coatran_eff_pop_size', required=True, type=float, help="CoaTran (Constant): Effective Population Size")
    parser.add_argument('--time_tree_seed_height', required=True, type=float, help="Time Tree: Seed Height (months)")
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
def load_sample_time_probs(sample_time_probs_fn, delim=','):
    probs = dict() # probs[(gender,risk,agerange,race,timetype)] = sampling probability
    header = True
    for l in open(sample_time_probs_fn):
        parts = [v.strip() for v in l.split(delim)]
        if header:
            header = False
            colnames = parts; name_to_ind = {colname:i for i,colname in enumerate(colnames)}
            for colname in SAMPLE_TIME_PROB_COLUMNS:
                if colname not in name_to_ind:
                    raise ValueError("Sample time probabilities (CSV) does not have a header called '%s': %s" % (colname,sample_time_probs_fn))
        else:
            k = tuple(parts[name_to_ind[colname]] for colname in SAMPLE_TIME_PROB_COLUMNS[:-1]) # [:-1] to skip 'probability'
            if k in probs:
                raise ValueError("Duplicate sample time probability encountered in sample time probabilities CSV: %s" % k)
            probs[k] = float(parts[name_to_ind['probability']])
    return probs

# load demographic data
def load_demographics(demographic_fn, delim='\t'):
    demographics = dict() # demographics[ID][category] = value
    header = True
    for l in open(demographic_fn):
        parts = [v.strip() for v in l.split(delim)]
        if header:
            header = False
            colnames = parts; name_to_ind = {colname:i for i,colname in enumerate(colnames)}
            for colname in DEMOGRAPHICS_COLUMNS:
                if colname not in name_to_ind:
                    raise ValueError("Demographic data output by abm_hiv-HRSA_SD does not have a column called '%s': %s" % (colname,demographic_fn))
        else:
            ID = parts[name_to_ind['id']]
            if ID in demographics:
                raise ValueError("Duplicate ID encountered in demographics file: %s" % ID)
            vals = {colname:parts[name_to_ind[colname]] for colname in DEMOGRAPHICS_COLUMNS[1:]} # [1:] to skip 'id'
            end_age = float(vals['age'])
            demographics[ID] = { # 'gender', 'risk', 'agerange', 'race'
                'gender': vals['gender'],
                'risk': vals['risk'],
                'race': vals['race'],
                'end_age': end_age,
            }
    return demographics

# get sample times from all possible times
def sample_times_from_all_times(outdir, end_time, demographic_fn, all_times_fn, probs, demographic_delim='\t', all_times_delim='\t'):
    demographics = load_demographics(demographic_fn, delim=demographic_delim)
    sample_times_fn = '%s/error_free_files/sample_times.txt' % outdir; f = open(sample_times_fn, 'w') # TSV file
    header = True
    for l in open(all_times_fn):
        if header:
            header = False
        else:
            ID, month, event = [v.strip() for v in l.split(all_times_delim)]; month = float(month)
            if ID not in demographics:
                raise KeyError("ID not found in demographic data: %s" % ID)
            demo = demographics[ID]
            k = [demo[colname] for colname in SAMPLE_TIME_PROB_COLUMNS[:-3]] # [:-2] to skip 'agerange', 'timetype', and 'probability
            age = int(demo['end_age'] + ((month - end_time) / 12)) # cast to int to round down
            agerange = None
            for ind, curr_range in enumerate(AGE_RANGES_FLOAT):
                if curr_range[0] <= age <= curr_range[1]:
                    agerange = AGE_RANGES[ind]; break
            if agerange is None:
                raise ValueError("Encountered age outside of valid age ranges: %s (ID: %s, month: %s, event: %s)" % (age,ID,month,event))
            k = tuple(k + [agerange, event])
            if k not in probs:
                raise KeyError("Probability not found: %s" % k)
            p = probs[k]
            if random() <= p:
                f.write("%s\t%s\n" % (ID,month))
    f.close()
    return sample_times_fn

def sample_time_tree(outdir, transmission_fn, sample_times_fn, eff_pop_size, seed_height, path_coatran_constant=DEFAULT_PATH_COATRAN_CONSTANT, verbose=True):
    command = [path_coatran_constant, transmission_fn, sample_times_fn, str(eff_pop_size)]
    if verbose:
        print_log("Running CoaTran (Constant) Command: %s" % ' '.join(command))
    coatran_stdout = check_output(command)
    time_trees_fn = '%s/error_free_files/phylogenetic_trees/separate_trees.time.tre' % outdir; f = open(time_trees_fn, 'wb'); f.write(coatran_stdout); f.close()
    if verbose:
        print_log("Separate time trees written to: %s" % time_trees_fn)
        print_log("Loading time trees into TreeSwift...")
    time_trees = [read_tree_newick(l) for l in coatran_stdout.decode().splitlines()]
    if verbose:
        print_log("Sampling seed time tree...")
    seed_time_tree = nonhomogeneous_yule_tree(lambda t: exp(-t**2)+1, end_num_leaves=len(time_trees))
    seed_time_tree.scale_edges(seed_height/seed_time_tree.height())
    seed_time_leaves = list(seed_time_tree.traverse_leaves())
    if verbose:
        print_log("Merging time trees into single time tree...")
    for i in range(max(len(time_trees), len(seed_time_leaves))): # max included to trigger index-out-of-bounds if somehow lengths are different
        seed_time_leaves[i].children.append(time_trees[i].root)
    if verbose:
        print_log("Suppressing unifurcations in merged time tree...")
    seed_time_tree.suppress_unifurcations()
    time_tree_fn = '%s/error_free_files/phylogenetic_trees/merged_tree.time.tre' % outdir
    seed_time_tree.write_tree_newick(time_tree_fn)
    return time_tree_fn

def scale_tree(outdir, time_tree_fn, mutation_rate_loc, mutation_rate_scale, verbose=True):
    if verbose:
        print_log("Reloading time tree...")
    tree = read_tree_newick(time_tree_fn)
    if verbose:
        print_log("Scaling time tree by sampling mutation rates...")
    nodes = list(tree.traverse_preorder())
    rates = truncnorm.rvs(a=0, b=float('inf'), loc=mutation_rate_loc, scale=mutation_rate_scale, size=len(nodes))
    for i in range(max(len(nodes), len(rates))):
        nodes[i].edge_length *= rates[i]
    mut_tree_fn = '%s/error_free_files/phylogenetic_trees/merged_tree.tre' % outdir
    tree.write_tree_newick(mut_tree_fn)
    return mut_tree_fn

# main execution
if __name__ == "__main__":
    # parse/check user args and set things up
    args = parse_args(); check_args(args)
    makedirs(args.output, exist_ok=True)
    makedirs('%s/error_free_files' % args.output, exist_ok=True)
    makedirs('%s/error_free_files/phylogenetic_trees' % args.output, exist_ok=True)
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
    print_log("Loading Sample Time Probabilities (CSV): %s" % args.sample_time_probs_csv)
    probs = load_sample_time_probs(args.sample_time_probs_csv)
    print_log()
    print_log("=== abm_hiv-HRSA_SD Progress ===")
    end_time, calibration_fn, transmission_fn, all_times_fn, demographic_fn = run_abm_hiv_hrsa_sd(
        args.output,                   # FAVITES-ABM-HIV-SD-Lite output directory
        args.abm_hiv_params_xlsx,      # parameter XLSX file
        args.abm_hiv_trans_start,      # starting transition rate
        args.abm_hiv_trans_end,        # ending transition rate
        args.abm_hiv_trans_time,       # time to go from starting to ending transition rate
        args.path_abm_hiv_commandline, # path to abm_hiv-HRSA_SD/abm_hiv_commandline.R script
        args.path_abm_hiv_modules)     # path to abm_hiv-HRSA_SD/modules folder
    print_log()
    print_log("Determining sample times...")
    sample_times_fn = sample_times_from_all_times(args.output, end_time, demographic_fn, all_times_fn, probs)
    print_log("Sample times written to: %s" % sample_times_fn)
    print_log(); print_log()

    # simulate time and mutation phylogenies
    print_log("===== PHYLOGENY =====")
    print_log("=== Time Tree Arguments ===")
    print_log("Effective Population Size: %s" % args.coatran_eff_pop_size)
    print_log("Seed Height (months): %s" % args.time_tree_seed_height)
    print_log()
    print_log("=== CoaTran (Constant) Progress ===")
    time_tree_fn = sample_time_tree(args.output, transmission_fn, sample_times_fn, args.coatran_eff_pop_size, args.time_tree_seed_height, DEFAULT_PATH_COATRAN_CONSTANT)
    print_log("Time tree written to: %s" % time_tree_fn)
    print_log()
    print_log("=== Mutation Tree Arguments ===")
    print_log("Mutation Rate Truncated Normal Location (mutations/month): %s" % args.mutation_rate_loc)
    print_log("Mutation Rate Truncated Normal Scale (mutations/month): %s" % args.mutation_rate_scale)
    print_log()
    print_log("=== Mutation Tree Progress ===")
    mut_tree_fn = scale_tree(args.output, time_tree_fn, args.mutation_rate_loc, args.mutation_rate_scale)
    print_log("Mutation tree written to: %s" % mut_tree_fn)
