#!/usr/bin/env python3
'''
Calculate various metrics from a mutation tree
'''

# imports
from csv import reader
from sys import argv
from treeswift import read_tree_newick
import networkx as nx

# useful constants
GENETIC_LINKAGE_THRESHOLD = 0.015

# load risk factors from demographics.csv
def load_risk_factors(demographics_csv_fn):
    # load all risk factors, but keep track of unique risk factors (for use in binary assortativity)
    all_risk_factors = dict(); unique_risk_factors = set()
    with open(demographics_csv_fn) as f:
        ind_ID = None; ind_risk = None
        for row_num, row in enumerate(reader(f)):
            if row_num == 0:
                for i, v in enumerate(row):
                    vup = v.strip().upper()
                    if vup in {'UCSD_ID', 'UCI'}:
                        ind_ID = i
                    elif vup in {'RISK', 'RISK_MT'}:
                        ind_risk = i
            else:
                curr_risk_factor = row[ind_risk].strip()
                all_risk_factors[row[ind_ID].strip()] = curr_risk_factor
                unique_risk_factors.add(curr_risk_factor)

    # create binary "YES/NO" risk factors and merge to return
    risk_factors = {'all': all_risk_factors}
    for curr_risk_factor in unique_risk_factors:
        risk_factors[curr_risk_factor] = {person:(curr_risk_factor == all_risk_factors[person]) for person in all_risk_factors}
    return risk_factors

# build genetic linkage network from distance matrix
def build_genetic_network(distances, risk_factors, start_year=0, only_new=False):
    # build genetic linkage network
    genetic_network = nx.Graph()
    for i in range(len(leaf_labels)-1):
        ul = leaf_labels[i]; u = ul.split('_')[0].strip(); uy = int(ul.split('_')[1].strip())
        for j in range(i+1, len(leaf_labels)):
            vl = leaf_labels[j]; v = vl.split('_')[0].strip(); vy = int(vl.split('_')[1].strip()); d = distances[ul][vl]
            if d <= GENETIC_LINKAGE_THRESHOLD and ((not only_new) or (uy >= start_year or vy >= start_year)):
                genetic_network.add_edge(u, v, weight=d)

    # label nodes with risk factors and return
    for curr_risk_factor in risk_factors:
        nx.set_node_attributes(genetic_network, risk_factors[curr_risk_factor], curr_risk_factor)
    return genetic_network

# calculate link proportions
def calc_link_proportions(genetic_network, risk_factors):
    risk_factor_labels = set(risk_factors.keys()) - {'all'}
    link_counts = {u:{v:0 for v in risk_factor_labels} for u in risk_factor_labels}
    for u, v in genetic_network.edges:
        risk_u = [l for l in risk_factor_labels if risk_factors[l][u]][0]
        risk_v = [l for l in risk_factor_labels if risk_factors[l][v]][0]
        link_counts[risk_u][risk_v] += 1
        link_counts[risk_v][risk_u] += 1
    return link_counts

# main execution
if __name__ == "__main__":
    # parse args and load data
    if len(argv) != 4:
        print("USAGE: %s <newick_mutation_tree> <demographics_csv> <start_year>" % argv[0]); exit(1)
    tree = read_tree_newick(argv[1])
    risk_factors = load_risk_factors(argv[2])
    start_year = int(argv[3])
    leaf_labels = [node.label for node in tree.traverse_leaves()]

    # generate molecular linkage network
    distances = tree.distance_matrix(leaf_labels=True)
    genetic_network = build_genetic_network(distances, risk_factors, only_new=False)
    link_proportions = calc_link_proportions(genetic_network, risk_factors)
    genetic_network_new = build_genetic_network(distances, risk_factors, start_year=start_year, only_new=True)
    link_proportions_new = calc_link_proportions(genetic_network_new, risk_factors)
    for u in link_proportions:
        for v in link_proportions[u]:
            print("Link Proportion (%s-%s): %s" % (u, v, link_proportions[u][v]/sum(link_proportions[u].values())))
    for u in link_proportions_new:
        for v in link_proportions_new[u]:
            print("New Link Proportion (%s-%s): %s" % (u, v, link_proportions_new[u][v]/sum(link_proportions_new[u].values())))
