#!/usr/bin/env python3
'''
Given a time tree and identical-topology mutation tree (Newick), calculate root-to-tip mutation rates
'''

# imports
from sys import argv
from treeswift import read_tree_newick, read_tree_nexus

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

# main execution
if __name__ == "__main__":
    if len(argv) != 3:
        print("USAGE: %s <time_tree> <mutation_tree>" % argv[0]); exit(1)
    time_tree = load_tree(argv[1])
    mut_tree = load_tree(argv[2])
    time_root_dists = {node.label.strip():d for node,d in time_tree.distances_from_root(leaves=True, internal=False)}
    mut_root_dists = {node.label.strip():d for node,d in mut_tree.distances_from_root(leaves=True, internal=False)}
    for l in time_root_dists:
        if l in mut_root_dists:
            print(mut_root_dists[l] / time_root_dists[l])
