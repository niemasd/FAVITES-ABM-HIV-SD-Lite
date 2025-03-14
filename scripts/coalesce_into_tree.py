#!/usr/bin/env python3
'''
Coalesce new sample times into an existing tree.

The coalescent 'rate' is the Poisson rate of a single pairwise coalescent event.
In other words, with N lineages, the overall coalescent rate is r * N(N-1)/2
'''
from numpy.random import exponential
from random import choice, random
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
    def discard(self, x):
        if x in self.indices:
            self.remove(x)
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
    if tree.root.edge_length is not None and tree.root.edge_length > 0:
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

    # sort nodes in reverse order of time and coalesce
    nodes_rev_time = sorted([(node_times[node], node) for node in node_times], reverse=True)
    curr_tree_nodes = RandomPopList(); curr_floating_nodes = RandomPopList()
    nodes_ind = 0; curr_time = nodes_rev_time[0][0]
    while nodes_ind < len(nodes_rev_time):
        # add next node encountered going backwards in time
        curr_time, next_node = nodes_rev_time[nodes_ind]; nodes_ind += 1
        if next_node.is_root(): # root of input tree, or floating new node
            curr_floating_nodes.insert(next_node)
        else:
            curr_tree_nodes.insert(next_node)
        for child in next_node.children:
            curr_tree_nodes.discard(child)

        # coalesce as much as possible
        while len(curr_floating_nodes) + len(curr_tree_nodes) != 1:
            #print('F:%s, T:%s, I:%s, L:%s' % (len(curr_floating_nodes), len(curr_tree_nodes), nodes_ind, len(nodes_rev_time)))
            # if no floating nodes right now, just skip to next time window
            if len(curr_floating_nodes) == 0:
                break

            # determine coalescent time (one node must be floating, the other can be floating or tree)
            num_pairs_ff = len(curr_floating_nodes) * (len(curr_floating_nodes)-1) / 2
            num_pairs_ft = len(curr_floating_nodes) * len(curr_tree_nodes)
            num_pairs_total = num_pairs_ff + num_pairs_ft
            coalescent_time = curr_time - exponential(scale=1./(rate*num_pairs_total))
            if nodes_ind < len(nodes_rev_time) and coalescent_time < nodes_rev_time[nodes_ind][0]:
                break

            # create new parent node and pick random floating node to coalesce
            p = Node(); node_times[p] = coalescent_time
            u = curr_floating_nodes.pop_random()
            p.add_child(u); u.edge_length = node_times[u] - node_times[p]

            # flip coin to see if coalescing with floating or tree node, and then randomly pick and coalesce
            if random() < len(curr_tree_nodes)/(len(curr_tree_nodes)+len(curr_floating_nodes)):
                v = curr_tree_nodes.pop_random(); curr_tree_nodes.insert(p)
                v.parent.add_child(p); v.parent.remove_child(v); p.edge_length = node_times[p] - node_times[p.parent]
            else:
                v = curr_floating_nodes.pop_random(); curr_floating_nodes.insert(p)
            p.add_child(v); v.edge_length = node_times[v] - node_times[p]
            curr_time = coalescent_time
            if u == tree.root or v == tree.root:
                tree.root = p

    # output tree
    print(tree.newick())
