#!/usr/bin/env python3
from csv import reader
from networkx import Graph
from os.path import isdir, isfile
from sys import argv
from treeswift import read_tree_newick
import re
GENETIC_LINKAGE_THRESHOLD = 0.015

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
def build_genetic_network(mutation_tree):
    distances = mutation_tree.distance_matrix(leaf_labels=True)
    leaf_labels = [node.label for node in mutation_tree.traverse_leaves()]
    genetic_network = Graph()
    for i in range(len(leaf_labels)-1):
        ul = leaf_labels[i]; u = ul.split('_')[0].strip()
        for j in range(i+1, len(leaf_labels)):
            vl = leaf_labels[j]; v = vl.split('_')[0].strip(); d = distances[ul][vl]
            if d <= GENETIC_LINKAGE_THRESHOLD:
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
def score(sim_out_folder, calibration_csv, verbose=True):
    calibration_data = load_calibration_data(calibration_csv)
    abm_calibration_data = load_abm_calibration_output('%s/abm_hiv_calibration_data.tsv' % sim_out_folder)
    abm_risk_factors = load_abm_demographics_output('%s/abm_hiv_demographic_data.tsv' % sim_out_folder)
    mutation_tree = read_tree_newick('%s/error_free_files/phylogenetic_trees/merged_tree.tre' % sim_out_folder)
    genetic_network = build_genetic_network(mutation_tree)
    link_proportions = calc_link_proportions(genetic_network, abm_risk_factors)
    score = 0
    sim_start_time = int(sorted(k for k in calibration_data.keys() if re.match(r'[0-9]{4}_',k))[0].split('_')[0])
    print("Calibration Key\tReal Value\tSimulation Value")
    for cal_key, cal_tup in calibration_data.items():
        sim_val = None; cal_val, cal_w = cal_tup
        if re.match(r'[0-9]{4}_', cal_key):
            try:
                year, risk = cal_key.split('_')
                year = int(year)
                risk = risk.replace(' & ', 'and')
                if risk.lower() == 'other':
                    risk = 'other' # 'other' is lowercase in the ABM R code output
                sim_val = abm_calibration_data['newinfects_agg'][year-sim_start_time][risk]
            except:
                pass # sim_val will remain None, so ValueError below will be thrown
        elif cal_key.startswith('Link_'):
            risk_u, risk_v = [x.strip().upper() for x in cal_key.split('_')[1].split('-')]
            sim_val = link_proportions[risk_u][risk_v] / sum(link_proportions[risk_u].values())
        if sim_val is None:
            raise ValueError("Unknown calibration key: %s" % cal_key)
        else:
            if verbose:
                print("%s\t%s\t%s" % (cal_key, cal_val, sim_val))
            score += cal_w * ((sim_val-cal_val)**2)
    score = score**0.5
    if verbose:
        print("Overall Calibration Score\tN/A\t%s" % score)
    return score

# main execution
if __name__ == "__main__":
    if len(argv) != 3:
        print("%s <sim_out_folder> <calibration_csv>" % argv[0]); exit(1)
    if not isdir(argv[1]):
        print("Folder not found: %s" % argv[1]); exit(1)
    if not isfile(argv[2]):
        print("File not found: %s" % argv[2]); exit(1)
    s = score(argv[1], argv[2])
