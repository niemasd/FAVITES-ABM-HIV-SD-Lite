#!/usr/bin/env python3
from csv import reader
from networkx import Graph
from os import remove
from os.path import isdir, isfile
from treeswift import read_tree_newick
import argparse
import re
GENETIC_LINKAGE_THRESHOLD = 0.015
CALIBRATION_MODES = {'epi', 'epi+genetic'}

# parse user args
def parse_args():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-s', '--simulation_folder', required=True, type=str, help="Simulation Output Folder")
    parser.add_argument('-c', '--calibration_csv', required=True, type=str, help="Calibration Spreadsheet (CSV)")
    parser.add_argument('-m', '--calibration_mode', required=True, type=str, help="Calibration Mode (%s)" % ', '.join(sorted(CALIBRATION_MODES)))
    parser.add_argument('-o', '--output', required=False, type=str, default='stdout', help="Scoring Output")
    args = parser.parse_args()
    if not isdir(args.simulation_folder):
        print("Folder not found: %s" % args.simulation_folder); exit(1)
    if not isfile(args.calibration_csv):
        print("File not found: %s" % args.calibration_csv); exit(1)
    args.calibration_mode = args.calibration_mode.strip().lower()
    if args.calibration_mode not in CALIBRATION_MODES:
        print("Invalid calibration mode: %s" % args.calibration_mode); exit(1)
    if isfile(args.output):
        print("Output file exists: %s" % args.output); exit(1)
    return args

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

# load demographic data output by the ABM R code
def load_abm_demographics_output(abm_demographics_tsv):
    risk_factors = dict()
    for row in reader(open(abm_demographics_tsv), delimiter='\t'):
        ID, gender, risk, age, race = [x.strip() for x in row]
        if ID == 'id':
            continue # skip header row
        risk_factors[ID] = risk.upper()
    return risk_factors

# build genetic linkage network from distance matrix
def build_genetic_network(mutation_tree, only_new=False):
    distances = mutation_tree.distance_matrix(leaf_labels=True)
    leaf_labels = [node.label for node in mutation_tree.traverse_leaves()]
    genetic_network = Graph()
    for i in range(len(leaf_labels)-1):
        ul = leaf_labels[i]; u = ul.split('_')[0].strip()
        for j in range(i+1, len(leaf_labels)):
            vl = leaf_labels[j]; v = vl.split('_')[0].strip(); d = distances[ul][vl]
            if d <= GENETIC_LINKAGE_THRESHOLD and ((not only_new) or ('|' in ul or '|' in vl)):
                genetic_network.add_edge(u, v, weight=d)
    return genetic_network

# calculate link proportions
def calc_link_proportions(genetic_network, abm_risk_factors):
    risk_factor_labels = set(abm_risk_factors.values())
    link_counts = {u:{v:0 for v in risk_factor_labels} for u in risk_factor_labels}
    for u, v in genetic_network.edges:
        if '|' in u:
            risk_u = abm_risk_factors[u.split('|')[1]]
        else:
            try:
                risk_u = abm_risk_factors[u.split('_')[0]]
            except:
                continue
        if '|' in v:
            risk_v = abm_risk_factors[v.split('|')[1]]
        else:
            try:
                risk_v = abm_risk_factors[v.split('_')[0]]
            except:
                continue
        link_counts[risk_u][risk_v] += 1
        link_counts[risk_v][risk_u] += 1
    return link_counts

# scoring function we want to optimize
def score(sim_out_folder, calibration_csv, calibration_mode, out_fn, verbose=True):
    if out_fn.strip().lower() == 'stdout':
        from sys import stdout as out_f
    else:
        out_f = open(out_fn, 'w')
    calibration_mode_parts = {v.strip().lower() for v in calibration_mode.split('+')}
    calibration_data = load_calibration_data(calibration_csv)
    abm_calibration_data = load_abm_calibration_output('%s/abm_hiv_calibration_data.tsv' % sim_out_folder)
    abm_risk_factors = load_abm_demographics_output('%s/abm_hiv_demographic_data.tsv' % sim_out_folder)
    mutation_tree = read_tree_newick('%s/error_free_files/phylogenetic_trees/merged_tree.tre' % sim_out_folder)
    genetic_network = None; link_proportions = None
    genetic_network_new = None; link_proportions_new = None
    score = 0
    sim_start_time = int(sorted(k for k in calibration_data.keys() if re.match(r'[0-9]{4}_',k))[0].split('_')[0])
    out_f.write("Calibration Key\tReal Value\tSimulation Value\n")
    for cal_key, cal_tup in calibration_data.items():
        sim_val = None; cal_val, cal_w = cal_tup
        if 'epi' in calibration_mode_parts:
            if re.match(r'[0-9]{4}_', cal_key):
                try:
                    year, risk = cal_key.split('_')
                    year = int(year)
                    risk = risk.replace(' & ', 'and')
                    if risk.lower() == 'other':
                        risk = 'other' # 'other' is lowercase in the ABM R code output
                    sim_val = abm_calibration_data['newinfects_agg'][year-sim_start_time][risk]
                except:
                    raise ValueError("Unknown calibration key: %s" % cal_key)
        if 'genetic' in calibration_mode_parts:
            if cal_key.startswith('Link_'):
                if genetic_network is None:
                    genetic_network = build_genetic_network(mutation_tree, only_new=False)
                if link_proportions is None:
                    link_proportions = calc_link_proportions(genetic_network, abm_risk_factors)
                try:
                    risk_u, risk_v = [x.strip().upper() for x in cal_key.split('_')[1].split('-')]
                    denominator = sum(link_proportions[risk_u].values())
                    if denominator == 0:
                        sim_val = 0
                    else:
                        sim_val = link_proportions[risk_u][risk_v] / sum(link_proportions[risk_u].values())
                except:
                    raise ValueError("Unknown calibration key: %s" % cal_key)
            if cal_key.startswith('NewLink_'):
                if genetic_network_new is None:
                    genetic_network_new = build_genetic_network(mutation_tree, only_new=True)
                if link_proportions_new is None:
                    link_proportions_new = calc_link_proportions(genetic_network_new, abm_risk_factors)
                try:
                    risk_u, risk_v = [x.strip().upper() for x in cal_key.split('_')[1].split('-')]
                    denominator = sum(link_proportions_new[risk_u].values())
                    if denominator == 0:
                        sim_val = 0
                    else:
                        sim_val = link_proportions_new[risk_u][risk_v] / sum(link_proportions_new[risk_u].values())
                except:
                    raise ValueError("Unknown calibration key: %s" % cal_key)
        if sim_val is not None:
            if verbose:
                out_f.write("%s\t%s\t%s\n" % (cal_key, cal_val, sim_val))
            score += cal_w * abs(sim_val-cal_val)
    if verbose:
        out_f.write("Overall Calibration Score\tN/A\t%s\n" % score)
    return score

# main execution
if __name__ == "__main__":
    args = parse_args()
    try:
        s = score(args.simulation_folder, args.calibration_csv, args.calibration_mode, args.output)
    except Exception as e:
        if isfile(args.output):
            remove(args.output) # delete output file if crash
        raise e
