#!/usr/bin/env python3
'''
Coalesce new sample times into an existing tree.

The coalescent 'rate' is the Poisson rate of a single pairwise coalescent event.
In other words, with N lineages, the overall coalescent rate is r * N(N-1)/2
'''
from numpy.random import exponential
from random import choice
from sys import argv
from treeswift import Node, read_tree_newick

# helper list class to randomly (not arbitrarily) pop elements
class RandomPopList:
    def __init__(self):
        self.elements = list()
        self.indices = dict()
    def __len__(self):
        return len(self.elements)
    def insert(self, x):
        if x in self.indices:
            raise ValueError("duplicate element: %s" % x)
        self.indices[x] = len(self.elements); self.elements.append(x)
    def remove(self, x):
        ind = self.indices[x]; del self.indices[x]; self.indices[self.elements[-1]] = ind
        self.elements[ind], self.elements[-1] = self.elements[-1], self.elements[ind]
        self.elements.pop()
    def pop_random(self):
        if len(self.elements) == 0:
            raise IndexError("pop from empty list")
        x = choice(self.elements); self.remove(x)
        return x

# main execution
if __name__ == "__main__":
    # parse args
    if '-h' in argv or '--help' in argv or len(argv) != 5:
        print("USAGE: %s <newick_tree> <tMRCA> <new_sample_times_tsv> <coalescent_rate>" % argv[0]); exit(1)
    tree = read_tree_newick(argv[1])
    root = tree.root
    if root.edge_length is not None and root.edge_length > 0:
        raise ValueError("Root cannot have an edge length")
    tmrca = float(argv[2])
    node_times = {Node(u.strip()):float(v) for u,v in [l.strip().split('\t') for l in open(argv[3])]} # (time, Node) tuples
    rate = float(argv[4])

    # add times of tree nodes
    for node in tree.traverse_preorder():
        if node.is_root():
            node_times[node] = tmrca
        else:
            node_times[node] = node_times[node.parent] + node.edge_length

    # sort nodes in reverse order of time
    nodes_rev_time = sorted([(node_times[node], node) for node in node_times], reverse=True)

    # coalesce nodes as much as possible
    curr_tree_nodes = RandomPopList(); curr_floating_nodes = RandomPopList()
    nodes_ind = 0; time = float('inf')
    while True:
        time = nodes_rev_time[nodes_ind][0]
        while nodes_ind < len(nodes_rev_time) and nodes_rev_time[nodes_ind][0] >= time:
            next_node = nodes_rev_time[nodes_ind][1]; nodes_ind += 1
            if next_node.parent is not None and next_node is not root:
                curr_floating_nodes.insert(next_node)
            else:
                curr_tree_nodes.insert(next_node)
                for child in next_node.children:
                    curr_tree_nodes.remove(child)
        if nodes_ind == len(nodes_rev_time):
            break # final coalescent window
        while len(curr_floating_nodes) != 0:
            n = len(curr_tree_nodes) + len(curr_floating_nodes)
            coalescent_time = time - exponential(scale=(2*rate)/(n*(n-1)))
            if coalescent_time < nodes_rev_time[nodes_ind][0]:
                break
            u = curr_floating_nodes.pop_random()
            pass # TODO pick v by flipping coin to see if v is tree or floating, then randomly pick v from that group
            pass # TODO coalesce u and v into new parent p, and add p to curr_tree_nodes

    # coalesce all remaining nodes
    while len(curr_floating_nodes) != 0:
        pass # TODO
