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
    for row_num, row in enumerate(reader(open(calibration_csv))):
        if row_num == 0:
            continue # skip header row
        d, v, w = [x.strip() for x in row]
        if d in calibration_data:
            raise ValueError("Duplicate entry in calibration CSV: %s" % d)
        calibration_data[d] = (float(v), float(w))
    return calibration_data

# load calibration data output by the ABM R code; data[metric][year][subgroup] = count
def load_abm_calibration_output(abm_calibration_tsv):
    data = dict()
    for row in reader(open(abm_calibration_tsv), delimiter='\t'):
        metric, demographic, month, subgroup, stat = [x.strip() for x in row]
        if metric.lower() == 'metric':
            continue # skip header row
        if metric not in data:
            data[metric] = dict()
        if demographic not in data[metric]:
            data[metric][demographic] = dict()
        year = (int(month) - 1) // 12
        if year not in data[metric][demographic]:
            data[metric][demographic][year] = dict()
        if subgroup not in data[metric][demographic][year]:
            data[metric][demographic][year][subgroup] = 0
        data[metric][demographic][year][subgroup] += int(stat)
    return data

# load demographic data output by the ABM R code
def load_abm_demographics_output(abm_demographics_tsv):
    demographics = dict()
    header = True
    for row in reader(open(abm_demographics_tsv), delimiter='\t'):
        if header:
            header = False
            name_to_ind = {colname:i for i,colname in enumerate(row)}
        else:
            ID = row[name_to_ind['id']]
            if ID in demographics:
                raise ValueError("Duplicate ID encountered in demographics file: %s" % ID)
            demographics[ID] = {colname:row[name_to_ind[colname]] for colname in name_to_ind if colname != 'id'}
            if 'age' in demographics[ID]: # end age
                demographics[ID]['age'] = float(demographics[ID]['age'])
    return demographics

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
def calc_link_proportions(genetic_network, abm_demographics, stratify_demographic):
    demographic_subgroups = {dems[stratify_demographic] for dems in abm_demographics.values()}
    link_counts = {u:{v:0 for v in demographic_subgroups} for u in demographic_subgroups}
    for u, v in genetic_network.edges:
        if '|' in u:
            subgroup_u = abm_demographics[u.split('|')[1]][stratify_demographic]
        else:
            try:
                subgroup_u = abm_demographics[u.split('_')[0]][stratify_demographic]
            except:
                continue
        if '|' in v:
            subgroup_v = abm_demographics[v.split('|')[1]][stratify_demographic]
        else:
            try:
                subgroup_v = abm_demographics[v.split('_')[0]][stratify_demographic]
            except:
                continue
        link_counts[subgroup_u][subgroup_v] += 1
        link_counts[subgroup_v][subgroup_u] += 1
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
    abm_demographics = load_abm_demographics_output('%s/abm_hiv_demographic_data.tsv' % sim_out_folder)
    mutation_tree = read_tree_newick('%s/error_free_files/phylogenetic_trees/merged_tree.tre' % sim_out_folder)
    genetic_network = None; link_proportions = dict()
    genetic_network_new = None; link_proportions_new = dict()
    score = 0
    sim_start_time = int(sorted(k for k in calibration_data.keys() if re.match(r'[0-9]{4}_',k))[0].split('_')[0])
    out_f.write("Calibration Key\tReal Value\tSimulation Value\n")
    for cal_key, cal_tup in calibration_data.items():
        sim_val = None; cal_val, cal_w = cal_tup

        # epi-specific calibration metrics
        if 'epi' in calibration_mode_parts:
            # YYYY_* calibration metrics
            if re.match(r'[0-9]{4}_', cal_key):
                try:
                    year, demographic, subgroup = cal_key.split('_') # e.g. 2019_race_black
                    year = int(year)
                    sim_val = abm_calibration_data['newinfects_agg'][demographic][year-sim_start_time][subgroup]
                except:
                    raise ValueError("Unknown calibration key: %s" % cal_key)

        # genetic-specific calibration metrics
        if 'genetic' in calibration_mode_parts:
            # Link_* calibration metrics = genetic linkages
            if cal_key.startswith('Link_'):
                DUMMY, demographic, uv = cal_key.split('_'); subgroup_u, subgroup_v = [x.strip() for x in uv.split('-')]
                if genetic_network is None:
                    genetic_network = build_genetic_network(mutation_tree, only_new=False)
                if demographic not in link_proportions:
                    link_proportions[demographic] = calc_link_proportions(genetic_network, abm_demographics, demographic)
                try:
                    denominator = sum(link_proportions[demographic][subgroup_u].values())
                    if denominator == 0:
                        sim_val = 0
                    else:
                        sim_val = link_proportions[demographic][subgroup_u][subgroup_v] / sum(link_proportions[subgroup_u].values())
                except:
                    raise ValueError("Unknown calibration key: %s" % cal_key)

            # NewLink_* calibration metrics = new genetic linkages
            if cal_key.startswith('NewLink_'):
                DUMMY, demographic, uv = cal_key.split('_'); subgroup_u, subgroup_v = [x.strip() for x in uv.split('-')]
                if genetic_network_new is None:
                    genetic_network_new = build_genetic_network(mutation_tree, only_new=True)
                if demographic not in link_proportions_new:
                    link_proportions_new[demographic] = calc_link_proportions(genetic_network_new, abm_demographics, demographic)
                try:
                    denominator = sum(link_proportions_new[demographic][subgroup_u].values())
                    if denominator == 0:
                        sim_val = 0
                    else:
                        sim_val = link_proportions_new[demographic][subgroup_u][subgroup_v] / sum(link_proportions_new[demographic][subgroup_u].values())
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
