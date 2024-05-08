#!/usr/bin/env python3
'''
Given a time tree and identical-topology mutation tree (Newick), calculate root-to-tip mutation rates
'''

# imports
from sys import argv
from treeswift import read_tree_newick

# main execution
if __name__ == "__main__":
    if len(argv) != 3:
        print("USAGE: %s <newick_time_tree> <newick_mutation_tree>" % argv[0]); exit(1)
    time_tree = read_tree_newick(argv[1])
    mut_tree = read_tree_newick(argv[2])
    time_root_dists = {node.label.strip():d for node,d in time_tree.distances_from_root(leaves=True, internal=False)}
    mut_root_dists = {node.label.strip():d for node,d in mut_tree.distances_from_root(leaves=True, internal=False)}
    for l in time_root_dists:
        if l in mut_root_dists:
            print(mut_root_dists[l] / time_root_dists[l])
