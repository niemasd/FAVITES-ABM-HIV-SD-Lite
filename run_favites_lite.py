#!/usr/bin/env python3
from datetime import datetime
from math import exp
from niemads import RandomSet
from numpy.random import exponential, lognormal
from openpyxl import load_workbook, Workbook
from openpyxl.cell import WriteOnlyCell
from openpyxl.styles import Alignment
from os import chdir, getcwd, makedirs
from os.path import abspath, expanduser, isdir, isfile
from random import random
from scipy.stats import gaussian_kde, truncexpon
from shutil import copy
from statistics import mean
from subprocess import check_output
from sys import argv, stdout
from time import mktime
from treesap.NHPP import NHPP_first_interarrival_time
from treeswift import Node, read_tree_newick, read_tree_nexus
from warnings import catch_warnings, simplefilter
import argparse

# useful constants
VERSION = '0.0.2'
LOGFILE = None
SAMPLE_TIME_PROB_COLUMNS = ['gender', 'risk', 'race', 'agerange', 'timetype', 'probability']
AGE_RANGES = ['0-24', '25-54', '55-999']
AGE_RANGES_FLOAT = [(float(v.split('-')[0]), float(v.split('-')[1])) for v in AGE_RANGES]
NHPP_RATE_FUNC = lambda t: exp(-t**2)+1
NHPP_EXP_RTT = 5.010476661748763

# defaults
DEFAULT_PATH_ABM_HIV_COMMANDLINE = "/usr/local/bin/abm_hiv-HRSA_SD/abm_hiv_commandline.R"
DEFAULT_PATH_ABM_HIV_MODULES = "/usr/local/bin/abm_hiv-HRSA_SD/modules"
DEFAULT_PATH_COATRAN_CONSTANT = "coatran_constant"
DEFAULT_FN_ABM_HIV_CALIBRATION = "abm_hiv_calibration_data.tsv"
DEFAULT_FN_ABM_HIV_DEMOGRAPHICS = "abm_hiv_demographic_data.tsv"
DEFAULT_FN_ABM_HIV_ID_MAP = 'abm_hiv_id_map.tsv'
DEFAULT_FN_ABM_HIV_LOG = "log_abm_hiv.txt"
DEFAULT_FN_ABM_HIV_TIMES = "abm_hiv_times.tsv"
DEFAULT_FN_LOG = "log_favites.txt"

# check if user is just printing version
if '--version' in argv:
    print("FAVITES-ABM-HIV-SD-Lite version %s" % VERSION); exit()

# return the current time as a string
def get_time():
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")

# print to the log (None implies stdout only)
def print_log(s='', end='\n'):
    tmp = "[%s] %s" % (get_time(), s)
    print(tmp, end=end); stdout.flush()
    if LOGFILE is not None:
        print(tmp, file=LOGFILE, end=end); LOGFILE.flush()

# load tree from file
def load_tree(fn):
    if fn.split('.')[-1].strip().lower() in {'nex', 'nexus'}:
        tmp = read_tree_nexus(fn)
        for k in ['taxlabels', 'translate']:
            if k in tmp:
                del tmp[k]
        if len(tmp) != 1:
            print("Nexus file must have exactly 1 tree: %s" % fn); exit(1)
        return tmp[list(tmp.keys())[0]]
    else:
        return read_tree_newick(fn)

# run abm_hiv-HRSA_SD
def run_abm_hiv_hrsa_sd(outdir, abm_hiv_params_xlsx, abm_hiv_sd_demographics_csv, abm_hiv_trans_start, abm_hiv_trans_end, abm_hiv_trans_time, path_abm_hiv_commandline=DEFAULT_PATH_ABM_HIV_COMMANDLINE, path_abm_hiv_modules=DEFAULT_PATH_ABM_HIV_MODULES, verbose=True):
    # run the R script
    command = ['Rscript', path_abm_hiv_commandline, path_abm_hiv_modules, abm_hiv_params_xlsx, abm_hiv_sd_demographics_csv, str(abm_hiv_trans_start), str(abm_hiv_trans_end), str(abm_hiv_trans_time)]
    if verbose:
        print_log("Running abm_hiv-HRSA_SD Command: %s" % ' '.join(command))
    log_f = open('%s/%s' % (outdir,DEFAULT_FN_ABM_HIV_LOG), 'w')
    log_f.write("=== ABM STDERR ===\n"); log_f.flush()
    abm_out = check_output(command, stderr=log_f).decode(); log_f.flush()
    log_f.write("\n\n=== ABM STDOUT ===\n"); log_f.write(abm_out); log_f.close()

    # separate the output
    if verbose:
        print_log("Parsing abm_hiv-HRSA_SD output...")
    abm_out = abm_out.split('[1] "Initializing population..."')[1]
    id_map_data, abm_out = abm_out.split('[1] "Starting simulation..."')
    sim_progress_data, abm_out = abm_out.split('[1] "Printing output..."')
    abm_out = abm_out.split('[1] "Calibration metrics..."')[1]
    calibration_data, abm_out = abm_out.strip().split('[1] "Transmission tree..."')
    transmission_data, abm_out = abm_out.strip().split('[1] "Sequence sample times..."')
    times_data, demographic_data = abm_out.strip().split('[1] "PLWH demographics..."')

    # parse simulation duration
    sim_duration = float(sim_progress_data.split('"Month:')[-1].split('"')[0]) / 12. # scale times from months to years
    if verbose:
        print_log("Simulation duration (years): %s" % sim_duration)

    # parse the ID map data and write to file
    if verbose:
        print_log("Parsing ID map data...")
    id_map_data = [[v.strip() for v in l.strip().split()] for l in id_map_data.strip().splitlines()]
    id_map_fn = '%s/%s' % (outdir, DEFAULT_FN_ABM_HIV_ID_MAP); f = open(id_map_fn, 'w')
    for row in id_map_data:
        if row[1].upper() != 'NA':
            f.write('\t'.join(row) + '\n')
    f.close()
    if verbose:
        print_log("ID map data written to: %s" % id_map_fn)

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
        u, v, t = [x.strip() for x in l.strip().split()]
        t = float(t) / 12. # scale times from months to years
        f.write('%s\t%s\t%s\n' % (u,v,t))
    f.close()
    if verbose:
        print_log("Transmission network written to: %s" % transmission_fn)

    # parse the times data and write to file
    if verbose:
        print_log("Parsing all event times...")
    times_data = [[v.strip() for v in l.strip().split()] for l in times_data.strip().splitlines()]
    all_times_fn = '%s/%s' % (outdir, DEFAULT_FN_ABM_HIV_TIMES); f = open(all_times_fn, 'w')
    for row_num, row in enumerate(times_data):
        if row_num == 0:
            f.write('ID\ttime (year)\tevent\n')
        else:
            f.write('%s\t%s\t%s\n' % (row[0], float(row[1]) / 12., row[2])) # scale times from months to years
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
    return sim_duration, calibration_fn, transmission_fn, all_times_fn, demographic_fn, id_map_fn

# parse user args
def parse_args():
    # arg parser
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-o', '--output', required=True, type=str, help="Output Directory")
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
    parser.add_argument('--gzip_output', action='store_true', help="Gzip Compress Output Files")
    parser.add_argument('--path_abm_hiv_commandline', required=False, type=str, default=DEFAULT_PATH_ABM_HIV_COMMANDLINE, help="Path to abm_hiv-HRSA_SD/abm_hiv_commandline.R")
    parser.add_argument('--path_abm_hiv_modules', required=False, type=str, default=DEFAULT_PATH_ABM_HIV_MODULES, help="Path to abm_hiv-HRSA_SD/modules")
    parser.add_argument('--path_coatran_constant', required=False, type=str, default=DEFAULT_PATH_COATRAN_CONSTANT, help="Path to coatran_constant")
    parser.add_argument('--abm_xlsx_rng_seed', required=False, type=int, default=None, help="Override ABM XLSX RNG Seed")
    parser.add_argument('--abm_xlsx_idu_injections_per_month', required=False, type=float, default=None, help="Override ABM XLSX IDU Injections Per Month")
    parser.add_argument('--abm_xlsx_het_partners_per_person_per_month', required=False, type=float, default=None, help="Override ABM XLSX Heterosexual Partners Per Person Per Month")
    parser.add_argument('--abm_xlsx_het_assortativity_risk', required=False, type=float, default=None, help="Override ABM XLSX Heterosexual Assortativity (Risk)")
    parser.add_argument('--abm_xlsx_het_assortativity_hiv_status', required=False, type=float, default=None, help="Override ABM XLSX Heterosexual Assortativity (HIV Status)")
    parser.add_argument('--abm_xlsx_msm_partners_per_person_per_month', required=False, type=float, default=None, help="Override ABM XLSX MSM Partners Per Person Per Month")
    parser.add_argument('--abm_xlsx_msm_assortativity_risk', required=False, type=float, default=None, help="Override ABM XLSX MSM Assortativity (Risk)")
    parser.add_argument('--abm_xlsx_msm_assortativity_hiv_status', required=False, type=float, default=None, help="Override ABM XLSX MSM Assortativity (HIV Status)")
    parser.add_argument('--abm_xlsx_msm_hiv_monthly_transmission_prob_mult', required=False, type=float, default=None, help="Override ABM XLSX MSM HIV Monthly Transmission Probability Between Pairs Multiplier")
    parser.add_argument('--abm_xlsx_msmw_percentage_msm_msmw', required=False, type=float, default=None, help="Override ABM XLSX MSMW Percentage of MSM who are MSMW")
    parser.add_argument('--version', action='store_true', help="Display Version")
    args = parser.parse_args()

    # fix/update args and return
    args.output = abspath(expanduser(args.output))
    args.abm_hiv_params_xlsx = abspath(expanduser(args.abm_hiv_params_xlsx))
    args.abm_hiv_sd_demographics_csv = abspath(expanduser(args.abm_hiv_sd_demographics_csv))
    args.sample_time_probs_csv = abspath(expanduser(args.sample_time_probs_csv))
    args.time_tree_seed = abspath(expanduser(args.time_tree_seed))
    args.mutation_tree_seed = abspath(expanduser(args.mutation_tree_seed))
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

    # check abm_hiv-HRSA_SD demographics CSV
    if not isfile(args.abm_hiv_sd_demographics_csv):
        raise ValueError("abm_hiv-HRSA_SD demographics CSV file not found: %s" % args.abm_hiv_sd_demographics_csv)

    # check starting transition rate
    if args.abm_hiv_trans_start <= 0:
        raise ValueError("abm_hiv-HRSA_SD starting transition rate must be positive: %s" % args.abm_hiv_trans_start)

    # check ending transition rate
    if args.abm_hiv_trans_end <= 0:
        raise ValueError("abm_hiv-HRSA_SD ending transition rate must be positive: %s" % args.abm_hiv_trans_end)

    # check time to go from starting to ending transition rate
    if args.abm_hiv_trans_time <= 0:
        raise ValueError("abm_hiv-HRSA_SD time from starting to ending transition rate must be positive: %s" % args.abm_hiv_trans_time)

    # check sample time probabilities
    if not isfile(args.sample_time_probs_csv):
        raise ValueError("Sample time probabilities (CSV) not found: %s" % args.sample_time_probs_csv)

    # check path_abm_hiv_commandline
    if not isfile(args.path_abm_hiv_commandline):
        raise ValueError("abm_hiv-HRSA_SD/abm_hiv_commandline.R not found: %s" % args.path_abm_hiv_commandline)

    # check path_abm_hiv_modules
    if not isdir(args.path_abm_hiv_modules):
        raise ValueError("abm_hiv-HRSA_SD/modules not found: %s" % args.path_abm_hiv_modules)

    # check seed real time tree
    if not isfile(args.time_tree_seed):
        raise ValueError("Seed time tree file not found: %s" % args.time_tree_seed)

    # check seed real mutation tree
    if not isfile(args.mutation_tree_seed):
        raise ValueError("Seed mutation tree file not found: %s" % args.mutation_tree_seed)

    # check tMRCA
    if args.time_tree_tmrca < 0:
        raise ValueError("Time Tree tMRCA must be positive")

    # check simulation start time
    if args.sim_start_time < args.time_tree_tmrca:
        raise ValueError("Simulation start time must be larger than time tree tMRCA")

    # check ABM XLSX overrides
    if args.abm_xlsx_rng_seed is not None and args.abm_xlsx_rng_seed < 0:
        raise ValueError("ABM XLSX RNG Seed must be non-negative: %s" % args.abm_xlsx_rng_seed)
    if args.abm_xlsx_idu_injections_per_month is not None and args.abm_xlsx_idu_injections_per_month < 0:
        raise ValueError("ABM XLSX IDU Injections Per Month must be non-negative: %s" % args.abm_xlsx_idu_injections_per_month)
    if args.abm_xlsx_het_partners_per_person_per_month is not None and args.abm_xlsx_het_partners_per_person_per_month < 0:
        raise ValueError("ABM XLSX Heterosexual Partners Per Person Per Month must be non-negative: %s" % args.abm_xlsx_het_partners_per_person_per_month)
    if args.abm_xlsx_msm_partners_per_person_per_month is not None and args.abm_xlsx_msm_partners_per_person_per_month < 0:
        raise ValueError("ABM XLSX MSM Partners Per Person Per Month must be non-negative: %s" % args.abm_xlsx_msm_partners_per_person_per_month)
    if args.abm_xlsx_msm_hiv_monthly_transmission_prob_mult is not None and args.abm_xlsx_msm_hiv_monthly_transmission_prob_mult < 0:
        raise ValueError("ABM XLSX MSM HIV Monthly Transmission Probability Between Partners Multiplier must be non-negative: %s" % args.abm_xlsx_msm_hiv_monthly_transmission_prob_mult)
    if args.abm_xlsx_msmw_percentage_msm_msmw is not None and (args.abm_xlsx_msmw_percentage_msm_msmw < 0 or args.abm_xlsx_msmw_percentage_msm_msmw > 1):
        raise ValueError("ABM XLSX MSMW Percentage of MSM who are MSMW must be in the range [0, 1]: %s" % args.abm_xlsx_msmw_percentage_msm_msmw)

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
            for colname in ['id', 'age']:
                if colname not in name_to_ind:
                    raise ValueError("Demographic data output by abm_hiv-HRSA_SD does not have a column called '%s': %s" % (colname,demographic_fn))
        else:
            ID = parts[name_to_ind['id']]
            if ID in demographics:
                raise ValueError("Duplicate ID encountered in demographics file: %s" % ID)
            demographics[ID] = {colname:parts[name_to_ind[colname]] for colname in name_to_ind if colname != 'id'}
            if 'age' in demographics[ID]: # end age
                demographics[ID]['age'] = float(demographics[ID]['age'])
    return demographics

# get sample times from all possible times
def sample_times_from_all_times(outdir, sim_duration, demographic_fn, all_times_fn, probs, demographic_delim='\t', all_times_delim='\t'):
    demographics = load_demographics(demographic_fn, delim=demographic_delim)
    sample_times_fn = '%s/error_free_files/sample_times.txt' % outdir; f = open(sample_times_fn, 'w') # TSV file
    header = True; already_sampled = set()
    for l in open(all_times_fn):
        if header:
            header = False
        else:
            ID, t, event = [v.strip() for v in l.split(all_times_delim)]; t = float(t)
            if ID in already_sampled:
                continue
            if ID not in demographics:
                raise KeyError("ID not found in demographic data: %s" % ID)
            demo = demographics[ID]
            k = [demo[colname] for colname in SAMPLE_TIME_PROB_COLUMNS[:-3]] # [:-2] to skip 'agerange', 'timetype', and 'probability
            age = int(demo['age'] + t - sim_duration) # cast to int to round down
            agerange = None
            for ind, curr_range in enumerate(AGE_RANGES_FLOAT):
                if curr_range[0] <= age <= curr_range[1]:
                    agerange = AGE_RANGES[ind]; break
            if agerange is None:
                raise ValueError("Encountered age outside of valid age ranges: %s (ID: %s, time: %s, event: %s)" % (age,ID,t,event))
            k = tuple(k + [agerange, event])
            if k in probs:
                p = probs[k]
                if random() <= p:
                    f.write("%s\t%s\n" % (ID,t)); already_sampled.add(ID)
    f.close()
    return sample_times_fn

def sample_time_tree(outdir, transmission_fn, sample_times_fn, id_map_fn, eff_pop_size, time_tree_seed_fn, time_tree_tmrca, sim_start_time, merge_model='yule', only_include_mapped=False, path_coatran_constant=DEFAULT_PATH_COATRAN_CONSTANT, verbose=True):
    # map individuals to their seed
    person_to_seed = dict(); infector = dict()
    for l in open(transmission_fn):
        u, v, t = [v.strip() for v in l.split('\t')]
        if v in person_to_seed:
            raise ValueError("Multi-infection in transmission network: %s" % v)
        if u == 'None':
            infector[v] = None; person_to_seed[v] = v
        else:
            infector[v] = u; person_to_seed[v] = person_to_seed[u]

    # run CoaTran to get a phylogeny for each seed's transmission chain
    command = [path_coatran_constant, transmission_fn, sample_times_fn, str(eff_pop_size)]
    if verbose:
        print_log("Running CoaTran (Constant) Command: %s" % ' '.join(command))
    coatran_stdout = check_output(command)
    time_trees_fn = '%s/error_free_files/phylogenetic_trees/separate_trees.time.tre' % outdir; f = open(time_trees_fn, 'wb'); f.write(coatran_stdout); f.close()
    if verbose:
        print_log("Separate time trees written to: %s" % time_trees_fn)

    # load time trees
    if verbose:
        print_log("Loading time trees into TreeSwift...")
    time_tree_objects = [read_tree_newick(l) for l in coatran_stdout.decode().splitlines()]
    time_trees = dict()
    for tree in time_tree_objects:
        seed = None
        for leaf in tree.traverse_leaves():
            person = leaf.label.split('|')[1].strip()
            if person not in person_to_seed:
                raise ValueError("Person not in transmission network: %s" % person)
            seed = person_to_seed[person]
            if seed in time_trees:
                raise ValueError("Multiple phylogenies for seed: %s" % seed)
            time_trees[seed] = tree; break
    if verbose:
        print_log("Loading seed time tree: %s" % time_tree_seed_fn)
    seed_time_tree = load_tree(time_tree_seed_fn)

    # merge time trees
    if verbose:
        print_log("Merging time trees into single time tree...")
    id_sim_to_real = {l.split('\t')[0].strip() : l.split('\t')[1].strip() for i,l in enumerate(open(id_map_fn)) if i != 0}
    merged_time_tree = merge_trees(outdir, seed_time_tree, time_trees, id_sim_to_real, time_tree_tmrca, sim_start_time, model=merge_model, only_include_mapped=only_include_mapped, verbose=verbose)
    if verbose:
        print_log("Suppressing unifurcations in merged time tree...")
    merged_time_tree.suppress_unifurcations()
    time_tree_fn = '%s/error_free_files/phylogenetic_trees/merged_tree.time.tre' % outdir
    merged_time_tree.write_tree_newick(time_tree_fn)
    return time_tree_fn

def merge_trees(outdir, seed_time_tree, time_trees, id_sim_to_real, time_tree_tmrca, sim_start_time, model='yule', only_include_mapped=False, verbose=True):
    # check that model is valid
    model = model.lower()
    if model not in {'yule', 'nhpp'}:
        raise ValueError("Invalid phylogenetic model: %s" % model)

    # label each node with its time and root distance (as # edges) for convenience/simplicity
    for node in seed_time_tree.traverse_preorder():
        if node.is_root():
            node.time = time_tree_tmrca
        else:
            # if the node has a label ending with _YYYYMMDD, use that as its time
            try:
                yyyymmdd = node.label.split('_')[-1]
                curr_date = datetime.strptime(yyyymmdd, '%Y%m%d')
                assert 1900 <= curr_date.year <= 2100
                # logic from: https://stackoverflow.com/a/6451892/2134991
                dt_curr_year_start = datetime(year=curr_date.year, month=1, day=1)
                dt_next_year_start = datetime(year=curr_date.year+1, month=1, day=1)
                year_elapsed = mktime(curr_date.timetuple()) - mktime(dt_curr_year_start.timetuple())
                year_duration = mktime(dt_next_year_start.timetuple()) - mktime(dt_curr_year_start.timetuple())
                year_elapsed_fraction = year_elapsed / year_duration
                node.time = curr_date.year + year_elapsed_fraction

            # otherwise, infer the node's time from its parent's time + incident branch length
            except:
                node.time = node.parent.time
                if node.edge_length is not None:
                    node.time += node.edge_length

    # set things up
    id_real_to_sim = {v:k for k,v in id_sim_to_real.items()}
    if only_include_mapped:
        sampled_seed_leaf_labels = {node.label for node in seed_time_tree.traverse_leaves() if node.time < sim_start_time and node.label.split('_')[0].strip() in id_real_to_sim}
    else:
        sampled_seed_leaf_labels = {node.label for node in seed_time_tree.traverse_leaves() if node.time < sim_start_time}
    merged_time_tree = seed_time_tree.extract_tree_with(sampled_seed_leaf_labels) # currently just seed tree with unsampled seeds removed
    merged_time_tree.suppress_unifurcations()
    for node in merged_time_tree.traverse_preorder():
        for k in ['node_params', 'edge_params']:
            if hasattr(node, k):
                delattr(node, k)
        if len(node.children) != 0:
            node.label = None
    base_time_tree_fn = '%s/error_free_files/phylogenetic_trees/base_tree.time.tre' % outdir
    merged_time_tree.write_tree_newick(base_time_tree_fn)
    merged_time_tree_rtt = {node:dist for node, dist in merged_time_tree.distances_from_root(leaves=True, internal=False, unlabeled=True, weighted=True)}
    seed_to_seed_leaf = {id_real_to_sim[node.label.split('_')[0].strip()]:node for node in merged_time_tree.traverse_leaves() if node.label.split('_')[0].strip() in id_real_to_sim}
    seeds_with_leaves = set(seed_to_seed_leaf.keys())
    seeds_without_leaves = RandomSet(set(time_trees.keys()) - seeds_with_leaves)
    merges_to_perform = [(seeds_without_leaves, merged_time_tree.root)]

    # label each node with its time and root distance (as # edges) for convenience/simplicity
    for tree in [merged_time_tree] + list(time_trees.values()):
        for node in tree.traverse_preorder():
            if node.is_root():
                node.edge_length = None; node.root_dist = 0
                if node is merged_time_tree.root:
                    node.time = time_tree_tmrca
                else:
                    node.time = sim_start_time
            else:
                if node.edge_length is None:
                    node.edge_length = 0
                node.time = node.parent.time + node.edge_length
                node.root_dist = node.parent.root_dist + 1

    # learn Yule splitting rate: expected BL = 1/(2lambda) -> lambda = 1/(2 * expected BL) = number of branches / (2 * sum of BL)
    if model == 'yule':
        tmp_bls = [node.edge_length for node in merged_time_tree.traverse_preorder() if not node.is_root()]
        yule_scale = 2. * sum(tmp_bls) / len(tmp_bls) # exponential scale = 1 / rate

    # set up NHPP model
    elif model == 'nhpp':
        avg_rtt = mean(merged_time_tree_rtt.values()); nhpp_scale_ratio = NHPP_EXP_RTT / avg_rtt

    # find roots of subtrees where to add seeds who do have a leaf
    for curr_seed, curr_seed_leaf in seed_to_seed_leaf.items():
        if curr_seed not in time_trees:
            continue
        curr_seed_tree = time_trees[curr_seed]
        avg_time_to_diagnosis = 1. # TODO REPLACE WITH DEMOGRAPHIC-STRATIFIED AVERAGE TIME TO DIAGNOSIS
        time_to_diagnosis = truncexpon.rvs(scale=avg_time_to_diagnosis, b=merged_time_tree_rtt[curr_seed_leaf])
        subtree_root = curr_seed_leaf
        while (not subtree_root.is_root()) and time_to_diagnosis >= subtree_root.edge_length:
            time_to_diagnosis -= subtree_root.edge_length; subtree_root = subtree_root.parent
        if subtree_root == curr_seed_leaf:
            subtree_root = curr_seed_leaf.parent
        merges_to_perform.append((RandomSet([curr_seed]), subtree_root))

    # find where to add each seed transmission chain phylogeny
    splits_to_add = dict() # keys = node objects in seed tree; values = list of phylo roots to add to the branch incident to this node
    for seeds_to_add, root in merges_to_perform:
        time_horizons = sorted(root.traverse_preorder(), key=lambda x: (x.time, x.root_dist))
        live_lineages = RandomSet([root]); horizon_ind = -1
        while True:
            horizon_ind += 1
            if horizon_ind == len(time_horizons) - 1: # skip latest node, as there won't be any lineages anymore
                horizon_ind = 0; live_lineages = RandomSet([root])
            curr_node = time_horizons[horizon_ind]; curr_time = curr_node.time; live_lineages.remove(curr_node)
            for child in curr_node.children:
                live_lineages.add(child)
            while True:
                if model == 'yule':
                    time_delta = exponential(scale=yule_scale/len(live_lineages))
                elif model == 'nhpp':
                    curr_nhpp_rate_func = lambda t: NHPP_RATE_FUNC((t + curr_time - time_tree_tmrca) * nhpp_scale_ratio)
                    rv = NHPP_first_interarrival_time(a=0); rv.set_L(curr_nhpp_rate_func)
                    time_delta = rv.rvs(size=1)[0] / nhpp_scale_ratio
                curr_time += time_delta
                if curr_time > time_horizons[horizon_ind+1].time:
                    break
                split_lineage = live_lineages.top_random()
                if split_lineage not in splits_to_add:
                    splits_to_add[split_lineage] = list()
                curr_seed = seeds_to_add.pop_random()
                curr_seed_root = time_trees[curr_seed].root
                curr_seed_root.insert_time = curr_time
                splits_to_add[split_lineage].append(curr_seed_root)
                if len(seeds_to_add) == 0:
                    break
            if len(seeds_to_add) == 0:
                break

    # actually perform the merge
    for node, roots in splits_to_add.items():
        for root in sorted(roots, key=lambda x: (x.insert_time, x.root_dist)):
            parent = node.get_parent(); parent.remove_child(node)
            new_node = Node(); new_node.time = root.insert_time
            new_node.add_child(root); new_node.add_child(node); parent.add_child(new_node)

    # fix branch lengths based on node times, and clean up tree
    merged_time_tree.suppress_unifurcations()
    for node in merged_time_tree.traverse_preorder():
        if node.is_root():
            node.edge_length = None # shouldn't be necessary, but just in case
        else:
            node.edge_length = node.time - node.parent.time
        if not node.is_leaf(): # delete internal labels
            node.label = None
    return merged_time_tree

# sample mutation rates learned from the real time_RTT / mut_RTT distribution
def scale_tree(outdir, sim_time_tree_fn, real_time_tree_fn, real_mutation_tree_fn, verbose=True):
    if verbose:
        print_log("Reloading time tree...")
    sim_tree = read_tree_newick(sim_time_tree_fn)
    if verbose:
        print_log("Loading real time and mutation trees...")
    real_time_tree = load_tree(real_time_tree_fn)
    real_mut_tree = load_tree(real_mutation_tree_fn)
    if verbose:
        print_log("Calculating time and mutation root-to-tip distances...")
    time_root_dists = {node.label.strip():d for node,d in real_time_tree.distances_from_root(leaves=True, internal=False)}
    mut_root_dists = {node.label.strip():d for node,d in real_mut_tree.distances_from_root(leaves=True, internal=False)}
    rtt_mut_rates = [mut_root_dists[l] / time_root_dists[l] for l in set(time_root_dists.keys()) & set(mut_root_dists.keys())]
    if verbose:
        print_log("Learning mutation rate distribution from root-to-tip distances...")
    mut_rate_dist = gaussian_kde(rtt_mut_rates)
    if verbose:
        print_log("Scaling time tree by sampling mutation rates...")
    nodes = [node for node in sim_tree.traverse_preorder() if not node.is_root()]
    rates = mut_rate_dist.resample(size=len(nodes))[0]
    for i in range(max(len(nodes), len(rates))):
        nodes[i].edge_length *= rates[i]
    mut_tree_fn = '%s/error_free_files/phylogenetic_trees/merged_tree.tre' % outdir
    sim_tree.write_tree_newick(mut_tree_fn)
    return mut_tree_fn

# main execution
if __name__ == "__main__":
    # parse/check user args and set things up
    args = parse_args(); check_args(args)
    makedirs(args.output, exist_ok=True)
    makedirs('%s/inputs' % args.output, exist_ok=True)
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

    # get ready for simulating contact and transmission network
    print_log("===== CONTACT AND TRANSMISSION NETWORK =====")
    print_log("=== Contact and Transmission Network Arguments ===")
    print_log("Simulation Start Time (year): %s" % args.sim_start_time)
    print_log("abm_hiv-HRSA_SD: Parameter XLSX: %s" % args.abm_hiv_params_xlsx)
    print_log("abm_hiv-HRSA_SD: Demographics CSV: %s" % args.abm_hiv_sd_demographics_csv)
    print_log("abm_hiv-HRSA_SD: Starting Transition Rate: %s" % args.abm_hiv_trans_start)
    print_log("abm_hiv-HRSA_SD: Ending Transition Rate: %s" % args.abm_hiv_trans_end)
    print_log("abm_hiv-HRSA_SD: Time (months) to Go from Starting to Ending Transition Rate: %s" % args.abm_hiv_trans_time)
    print_log("Loading Sample Time Probabilities (CSV): %s" % args.sample_time_probs_csv)
    probs_csv_copy_fn = '%s/inputs/sample_time_probs.csv' % args.output
    copy(args.sample_time_probs_csv, probs_csv_copy_fn)
    probs = load_sample_time_probs(probs_csv_copy_fn)
    print_log()
    print_log("=== abm_hiv-HRSA_SD Progress ===")

    # update data.xlsx file (if needed)
    data_xlsx_copy_fn = '%s/inputs/data.xlsx' % args.output
    wb = load_workbook(args.abm_hiv_params_xlsx, data_only=True)
    if args.abm_xlsx_rng_seed is not None:
        print_log("Overriding ABM RNG seed: %s" % args.abm_xlsx_rng_seed)
        wb['High Level Pop + Sim Features']['B4'] = args.abm_xlsx_rng_seed
    if args.abm_xlsx_idu_injections_per_month is not None:
        print_log("Overriding ABM IDU Injections Per Month: %s" % args.abm_xlsx_idu_injections_per_month)
        wb['Transmission + At-Risk Features']['B7'] = args.abm_xlsx_idu_injections_per_month
    if args.abm_xlsx_het_partners_per_person_per_month is not None:
        print_log("Overriding ABM Heterosexual Partners Per Person Per Month: %s" % args.abm_xlsx_het_partners_per_person_per_month)
        wb['Transmission + At-Risk Features']['D13'] = args.abm_xlsx_het_partners_per_person_per_month
    if args.abm_xlsx_het_assortativity_risk is not None:
        print_log("Overriding ABM Heterosexual Assortativity (Risk): %s" % args.abm_xlsx_het_assortativity_risk)
        wb['Transmission + At-Risk Features']['D17'] = args.abm_xlsx_het_assortativity_risk
    if args.abm_xlsx_het_assortativity_hiv_status is not None:
        print_log("Overriding ABM Heterosexual Assortativity (HIV Status): %s" % args.abm_xlsx_het_assortativity_hiv_status)
        wb['Transmission + At-Risk Features']['D18'] = args.abm_xlsx_het_assortativity_hiv_status
    if args.abm_xlsx_msm_partners_per_person_per_month is not None:
        print_log("Overriding ABM MSM Partners Per Person Per Month: %s" % args.abm_xlsx_msm_partners_per_person_per_month)
        wb['Transmission + At-Risk Features']['D30'] = args.abm_xlsx_msm_partners_per_person_per_month
    if args.abm_xlsx_msm_assortativity_risk is not None:
        print_log("Overriding ABM MSM Assortativity (Risk): %s" % args.abm_xlsx_msm_assortativity_risk)
        wb['Transmission + At-Risk Features']['D34'] = args.abm_xlsx_msm_assortativity_risk
    if args.abm_xlsx_msm_assortativity_hiv_status is not None:
        print_log("Overriding ABM MSM Assortativity (HIV Status): %s" % args.abm_xlsx_msm_assortativity_hiv_status)
        wb['Transmission + At-Risk Features']['D35'] = args.abm_xlsx_msm_assortativity_hiv_status
    if args.abm_xlsx_msm_hiv_monthly_transmission_prob_mult is not None:
        print_log("Overriding ABM MSM HIV Monthly Transmission Probability Between Partners Multiplier: %s" % args.abm_xlsx_msm_hiv_monthly_transmission_prob_mult)
        wb['Transmission + At-Risk Features']['D37'] = args.abm_xlsx_msm_hiv_monthly_transmission_prob_mult
    if args.abm_xlsx_msmw_percentage_msm_msmw is not None:
        print_log("Overriding ABM MSMW Percentage of MSM who are MSMW: %s" % args.abm_xlsx_msmw_percentage_msm_msmw)
        wb['Transmission + At-Risk Features']['D46'] = args.abm_xlsx_msmw_percentage_msm_msmw

    # fix data.xlsx and save
    for ws in wb.worksheets:
        for row in ws:
            for cell in row:
                if isinstance(cell.value, str):
                    cell.value = cell.value.replace('\n', '\r\n')
    wb.save(data_xlsx_copy_fn); wb.close()

    # run ABM R code
    dem_csv_copy_fn = '%s/inputs/demographics.csv' % args.output
    copy(args.abm_hiv_sd_demographics_csv, dem_csv_copy_fn)
    sim_duration, calibration_fn, transmission_fn, all_times_fn, demographic_fn, id_map_fn = run_abm_hiv_hrsa_sd(
        args.output,                   # FAVITES-ABM-HIV-SD-Lite output directory
        data_xlsx_copy_fn,             # parameter XLSX file
        dem_csv_copy_fn,               # demographics CSV file
        args.abm_hiv_trans_start,      # starting transition rate
        args.abm_hiv_trans_end,        # ending transition rate
        args.abm_hiv_trans_time,       # time to go from starting to ending transition rate
        args.path_abm_hiv_commandline, # path to abm_hiv-HRSA_SD/abm_hiv_commandline.R script
        args.path_abm_hiv_modules)     # path to abm_hiv-HRSA_SD/modules folder
    print_log()
    print_log("Determining sample times...")
    sample_times_fn = sample_times_from_all_times(args.output, sim_duration, demographic_fn, all_times_fn, probs)
    print_log("Sample times written to: %s" % sample_times_fn)
    print_log(); print_log()

    # simulate time and mutation phylogenies
    print_log("===== PHYLOGENY =====")
    print_log("=== Time Tree Arguments ===")
    print_log("Effective Population Size: %s" % args.coatran_eff_pop_size)
    print_log("Seed Time Tree: %s" % args.time_tree_seed)
    print_log("Seed Time Tree tMRCA: %s" % args.time_tree_tmrca)
    print_log()
    print_log("=== Mutation Tree Arguments ===")
    print_log("Seed Mutation Tree: %s" % args.mutation_tree_seed)
    print_log()
    print_log("=== CoaTran (Constant) Progress ===")
    copy_fns = list()
    for input_tree_fn, copy_fn_prefix in [(args.time_tree_seed, 'time_tree_seed'), (args.mutation_tree_seed, 'mutation_tree_seed')]:
        if input_tree_fn.split('.')[-1].strip().lower() in {'nexus', 'nex'}:
            copy_fn_ext = 'nex'
        else:
            copy_fn_ext = 'nwk'
        copy_fn = '%s/inputs/%s.%s' % (args.output, copy_fn_prefix, copy_fn_ext)
        copy(input_tree_fn, copy_fn)
        copy_fns.append(copy_fn)
    time_tree_seed_copy_fn = copy_fns[0]
    with catch_warnings():
        simplefilter("ignore")
        time_tree_fn = sample_time_tree(args.output, transmission_fn, sample_times_fn, id_map_fn, args.coatran_eff_pop_size, time_tree_seed_copy_fn, args.time_tree_tmrca, args.sim_start_time, merge_model='yule', only_include_mapped=args.time_tree_only_include_mapped, path_coatran_constant=args.path_coatran_constant)
    print_log("Time tree written to: %s" % time_tree_fn)
    print_log()
    print_log("=== Mutation Tree Progress ===")
    mut_tree_fn = scale_tree(args.output, time_tree_fn, args.time_tree_seed, args.mutation_tree_seed)
    print_log("Mutation tree written to: %s" % mut_tree_fn)
