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
    risk_factors = dict()
    with open(demographics_csv_fn) as f:
        ind_ID = None; ind_risk = None
        for row_num, row in enumerate(reader(f)):
            if row_num == 0:
                for i, v in enumerate(row):
                    vup = v.strip().upper()
                    if vup == 'UCI':
                        ind_ID = i
                    elif vup == 'EXPOSURE CATEGORY':
                        ind_risk = i
            else:
                risk_factors[row[ind_ID].strip()] = row[ind_risk].strip()
    return risk_factors

# build genetic linkage network from distance matrix
def build_genetic_network(distances, risk_factors):
    # build genetic linkage network
    genetic_network = nx.Graph()
    for i in range(len(leaf_labels)-1):
        ul = leaf_labels[i]; u = ul.split('_')[0].strip()
        for j in range(i+1, len(leaf_labels)):
            vl = leaf_labels[j]; v = vl.split('_')[0].strip(); d = distances[ul][vl]
            if d <= GENETIC_LINKAGE_THRESHOLD:
                genetic_network.add_edge(u, v, weight=d)

    # label nodes with risk factors and return
    nx.set_node_attributes(genetic_network, risk_factors, "risk_factor")
    return genetic_network

# main execution
if __name__ == "__main__":
    # parse args and load data
    if len(argv) != 3:
        print("USAGE: %s <newick_mutation_tree> <demographics_csv>" % argv[0]); exit(1)
    tree = read_tree_newick(argv[1])
    risk_factors = load_risk_factors(argv[2])
    leaf_labels = [node.label for node in tree.traverse_leaves()]

    # generate molecular linkage network
    distances = tree.distance_matrix(leaf_labels=True)
    genetic_network = build_genetic_network(distances, risk_factors)
    print("Assortativity (Risk Factor): %s" % nx.attribute_assortativity_coefficient(genetic_network, "risk_factor"))
