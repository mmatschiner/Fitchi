#!/usr/local/bin/python3

# Michael Matschiner, 2015-10-15
# michaelmatschiner@mac.com

# Import libraries and make sure we're on python 3.
import sys
if sys.version_info[0] < 3:
    print('ERROR: Python 3 is needed to run this script!')
    sys.exit(1)
import argparse
import textwrap
import random
import tempfile
import os
import re
import math
import scipy
import pygraphviz
from scipy import special
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.Alphabet import generic_dna

######################### networkx 2.0 Graph class source below #########################

#    Copyright (C) 2004-2015 by
#    Aric Hagberg <hagberg@lanl.gov>
#    Dan Schult <dschult@colgate.edu>
#    Pieter Swart <swart@lanl.gov>
#    All rights reserved.
#    BSD license.

class Graph(object):
    node_dict_factory = dict
    adjlist_dict_factory = dict
    edge_attr_dict_factory = dict

    def __init__(self, data=None, **attr):
        self.node_dict_factory = ndf = self.node_dict_factory
        self.adjlist_dict_factory = self.adjlist_dict_factory
        self.edge_attr_dict_factory = self.edge_attr_dict_factory

        self.graph = {}   # dictionary for graph attributes
        self.node = ndf()  # empty node attribute dict
        self.adj = ndf()  # empty adjacency dict
        # attempt to load graph with data
        if data is not None:
            convert.to_networkx_graph(data, create_using=self)
        # load graph attributes (must be after convert)
        self.graph.update(attr)
        self.edge = self.adj

    @property
    def name(self):
        return self.graph.get('name', '')

    @name.setter
    def name(self, s):
        self.graph['name'] = s

    def __iter__(self):
        return iter(self.node)

    def add_node(self, n, attr_dict=None, **attr):
        # set up attribute dict
        if attr_dict is None:
            attr_dict = attr
        else:
            try:
                attr_dict.update(attr)
            except AttributeError:
                raise NetworkXError(
                    "The attr_dict argument must be a dictionary.")
        if n not in self.node:
            self.adj[n] = self.adjlist_dict_factory()
            self.node[n] = attr_dict
        else:  # update attr even if node already exists
            self.node[n].update(attr_dict)

    def nodes_iter(self, data=False):
        if data:
            return iter(self.node.items())
        return iter(self.node)

    def nodes(self, data=False):
        return list(self.nodes_iter(data=data))

    def add_edge(self, u, v, attr_dict=None, **attr):
        # set up attribute dictionary
        if attr_dict is None:
            attr_dict = attr
        else:
            try:
                attr_dict.update(attr)
            except AttributeError:
                raise NetworkXError(
                    "The attr_dict argument must be a dictionary.")
        # add nodes
        if u not in self.node:
            self.adj[u] = self.adjlist_dict_factory()
            self.node[u] = {}
        if v not in self.node:
            self.adj[v] = self.adjlist_dict_factory()
            self.node[v] = {}
        # add the edge
        datadict = self.adj[u].get(v, self.edge_attr_dict_factory())
        datadict.update(attr_dict)
        self.adj[u][v] = datadict
        self.adj[v][u] = datadict

    def edges_iter(self, nbunch=None, data=False, default=None):
        seen = {}     # helper dict to keep track of multiply stored edges
        if nbunch is None:
            nodes_nbrs = self.adj.items()
        else:
            nodes_nbrs = ((n, self.adj[n]) for n in self.nbunch_iter(nbunch))
        if data is True:
            for n, nbrs in nodes_nbrs:
                for nbr, ddict in nbrs.items():
                    if nbr not in seen:
                        yield (n, nbr, ddict)
                seen[n] = 1
        elif data is not False:
            for n, nbrs in nodes_nbrs:
                for nbr, ddict in nbrs.items():
                    if nbr not in seen:
                        d = ddict[data] if data in ddict else default
                        yield (n, nbr, d)
                seen[n] = 1
        else:  # data is False
            for n, nbrs in nodes_nbrs:
                for nbr in nbrs:
                    if nbr not in seen:
                        yield (n, nbr)
                seen[n] = 1
        del seen

######################### networkx 2.0 Graph class source above #########################


# The Tree class.
class Tree(object):

    def __init__(self, newick_string):
        self.newick_string = newick_string
        self.nodes = []
        self.edges = []
        self.max_dist_to_root = 0
        self.max_number_of_edges = 0
        self.is_positioned = False
        self.pops = []

    def get_newick_string(self):
        return self.newick_string

    def parse_newick_string(self, pops):
        if pops == None:
            self.pops = []
        else:
            self.pops = pops
        number_of_internal_nodes = 0
        number_of_edges = 0

        # Remove comments from the tree string.
        pattern = re.compile("\[.*?\]")
        hit = "placeholder"
        while hit != None:
            hit = pattern.search(self.newick_string)
            if hit != None:
                self.newick_string = self.newick_string.replace(hit.group(0),"")

        # Check whether a branch above the root is present, and if so, remove it.
        if self.newick_string[0:2] == "((" and self.newick_string[-1] == ")" and self.newick_string[-2] != ")":
            level = 0
            newick_string_tail_start_pos = 0
            newick_string_tail = ""
            for pos in range(len(self.newick_string)):
                if self.newick_string[pos] == "(":
                    level += 1
                if self.newick_string[pos] == ")":
                    level -= 1
                if level == 1 and pos > 1:
                    newick_string_tail_start_pos = pos
                    newick_string_tail = self.newick_string[pos+1:]
                    break
            if newick_string_tail.count(",") == 0:
                self.newick_string = self.newick_string[1:newick_string_tail_start_pos+1]

        # Parse the bifurcating part of the tree.
        if ":" in self.newick_string:
            pattern = re.compile("\(([a-zA-Z0-9_\.\-]+?):[\d\.Ee-]+?,([a-zA-Z0-9_\.\-]+?):[\d\.Ee-]+?\)")
        else:
            pattern = re.compile("\(([a-zA-Z0-9_\.\-]+?),([a-zA-Z0-9_\.\-]+?)\)")
        hit = "placeholder"
        while hit != None:
            hit = pattern.search(self.newick_string)
            if hit != None:
                number_of_internal_nodes += 1
                number_of_edges += 2
                internal_node_id = "internalNode" + str(number_of_internal_nodes) + "X"
                edge1 = Edge("edge" + str(number_of_edges-1) + "X")
                self.increase_max_number_of_edges()
                edge1.set_node_ids([internal_node_id, hit.group(1)])
                node1 = Node(hit.group(1), False, self.pops)
                node1.set_parent_id(internal_node_id)
                if hit.group(1)[:12] == 'internalNode':
                    node1.set_size(0)
                else:
                    node1.set_size(1)
                    node1.add_record_id(hit.group(1))
                edge2 = Edge("edge" + str(number_of_edges) + "X")
                self.increase_max_number_of_edges()
                edge2.set_node_ids([internal_node_id, hit.group(2)])
                node2 = Node(hit.group(2), False, self.pops)
                node2.set_parent_id(internal_node_id)
                if hit.group(2)[:12] == 'internalNode':
                    node2.set_size(0)
                else:
                    node2.set_size(1)
                    node2.add_record_id(hit.group(2))
                self.edges.append(edge1)
                self.edges.append(edge2)
                self.nodes.append(node1)
                self.nodes.append(node2)
                self.newick_string = self.newick_string.replace(hit.group(0), internal_node_id)

        # Parse the remaining string with three branches (if the tree is unrooted).
        if ":" in self.newick_string:
            pattern_unrooted = re.compile("^\(([a-zA-Z0-9_\.\-]+?):[\d\.Ee-]+?,([a-zA-Z0-9_\.\-]+?):[\d\.Ee-]+?,([a-zA-Z0-9_\.\-]+?):[\d\.Ee-]+?\)$")
        else:
            pattern_unrooted = re.compile("^\(([a-zA-Z0-9_\.\-]+?),([a-zA-Z0-9_\.\-]+?),([a-zA-Z0-9_\.\-]+?)\)$")
        hit_unrooted = pattern_unrooted.search(self.newick_string)
        pattern_rooted = re.compile("^internalNode\d+X$")
        hit_rooted = pattern_rooted.search(self.newick_string)
        if hit_unrooted != None:
            number_of_internal_nodes += 1
            number_of_edges += 3
            root_node_id = "internalNode" + str(number_of_internal_nodes) + "X"
            edge1 = Edge("edge" + str(number_of_edges-2) + "X")
            self.increase_max_number_of_edges()
            edge1.set_node_ids([root_node_id, hit_unrooted.group(1)])
            edge2 = Edge("edge" + str(number_of_edges-1) + "X")
            self.increase_max_number_of_edges()
            edge2.set_node_ids([root_node_id, hit_unrooted.group(2)])
            edge3 = Edge("edge" + str(number_of_edges) + "X")
            self.increase_max_number_of_edges()
            edge3.set_node_ids([root_node_id, hit_unrooted.group(3)])
            node1 = Node(hit_unrooted.group(1), False, self.pops)
            node1.set_parent_id(root_node_id)
            node2 = Node(hit_unrooted.group(2), False, self.pops)
            node2.set_parent_id(root_node_id)
            node3 = Node(hit_unrooted.group(3), False, self.pops)
            node3.set_parent_id(root_node_id)
            node4 = Node(root_node_id, True, self.pops)
            node4.set_parent_id('None')
            if hit_unrooted.group(1)[:12] == 'internalNode':
                node1.set_size(0)
            else:
                node1.set_size(1)
            if hit_unrooted.group(2)[:12] == 'internalNode':
                node2.set_size(0)
            else:
                node2.set_size(1)
            if hit_unrooted.group(3)[:12] == 'internalNode':
                node3.set_size(0)
            else:
                node3.set_size(1)
            node4.set_size(0)
            self.edges.append(edge1)
            self.edges.append(edge2)
            self.edges.append(edge3)
            self.nodes.append(node1)
            self.nodes.append(node2)
            self.nodes.append(node3)
            self.nodes.append(node4)
        elif hit_rooted != None:
            node = Node(hit_rooted.group(0), True, self.pops)
            node.set_parent_id('None')
            node.set_size(0)
            self.nodes.append(node)
        else:
            print('ERROR: The newick tree string could not be parsed!')
            print(self.newick_string)
            sys.exit(1)

        # for node in self.nodes:
        #     print(node.info())

        # Add info about children to each node.
        count = 0
        for node in self.nodes:
            parent_found = False
            for parent in self.nodes:
                if node.get_parent_id() == parent.get_id():
                    parent.add_child_id(node.get_id())
                    parent_found = True
            if node.get_is_root() == False and parent_found == False:
                print("ERROR: The parent of node " + node.id + " named " + node.get_parent_id() + " could not be found!")
                sys.exit(1)

        # Calculate the distances to the root for each node.
        self.set_node_distances_to_root()

    def set_node_distances_to_root(self):
        # Reset max_dist_to_root.
        self.max_dist_to_root = 0
        # First make sure that one of the nodes is the root.
        there_is_a_root = False
        for node in self.nodes:
            if node.get_is_root():
                there_is_a_root = True
        if there_is_a_root == False:
            print("ERROR: No node is the root!")
            sys.exit(1)
        # Add distance to root for each node.
        for node in self.nodes:
            dist = 0
            if node.get_is_root() == False:
                root_found = False
                tmp_node = node
                while root_found == False:
                    dist += 1
                    # Find the parental node.
                    parent_found = False
                    for parent in self.nodes:
                        if tmp_node.get_parent_id() == parent.get_id():
                            tmp_node = parent
                            parent_found = True
                            break
                    if tmp_node.get_is_root():
                        root_found = True
                    elif parent_found == False:
                        print("ERROR: The parent of node " + tmp_node.id + " could not be found!")
                        sys.exit(1)

            node.set_distance_to_root(dist)
            if dist > self.max_dist_to_root:
                self.max_dist_to_root = dist

    def get_nodes(self):
        return self.nodes

    def get_edges(self):
        return self.edges

    def get_number_of_nodes(self):
        return len(self.nodes)

    def get_number_of_edges(self):
        return len(self.edges)

    def get_max_number_of_edges(self):
        return self.max_number_of_edges

    def increase_max_number_of_edges(self):
        self.max_number_of_edges += 1

    def get_max_dist_to_root(self):
        return self.max_dist_to_root

    def reconstruct_ancestral_sequences(self):
        # This uses the Fitch algorithm to assign sequences to internal nodes.

        # Find the sequence length.
        seq_length = 0
        for node in self.nodes:
            sequences = node.get_sequences()
            if sequences != []:
                sequence = sequences[0]
                seq_length = len(sequence)
        if seq_length == 0:
            print("WARNING: No sequence was found for any of the nodes!")

        # Prepare empty state sets and states of all nodes.
        for node in self.nodes:
            node.prepare_state_sets(seq_length)
            node.prepare_states(seq_length)

        # Compute the state sets of all nodes, sorted by their distance to the root (bottom up).
        invest_dist = self.max_dist_to_root
        while invest_dist >= 0:
            for node in self.nodes:
                if node.get_distance_to_root() == invest_dist:
                    node_id = node.get_id()
                    # If the node is an internal node, find its children and copy their state sets.
                    if node_id[:12] == 'internalNode':
                        children = []
                        for child in self.nodes:
                            if child.get_parent_id() == node_id:
                                children.append(child)
                        for x in range(seq_length):
                            children_state_sets = []
                            for child in children:
                                children_state_sets.append(child.get_state_set(x))
                            if len(children_state_sets) == 2:
                                intersection = list(set.intersection(set(children_state_sets[0]),set(children_state_sets[1])))
                            elif len(children_state_sets) == 3:
                                intersection = list(set.intersection(set(children_state_sets[0]),set(children_state_sets[1]),set(children_state_sets[2])))
                            else:
                                print('WARNING: An unexpected number of child state sets was encountered.')
                            if len(intersection) == 0:
                                # Flatten children state sets
                                flat_children_state_sets = []
                                for children_state_set in children_state_sets:
                                    for children_state in children_state_set:
                                        flat_children_state_sets.append(children_state)
                                state_set = list(set(flat_children_state_sets))
                            else:
                                state_set = intersection
                            node.set_state_set(x, state_set)

                    # If the node is terminal, get its state sets directly from its sequence.
                    else:
                        for x in range(seq_length):
                            if node.get_sequences()[0][x].upper() == 'A':
                                state_set = ['A']
                            elif node.get_sequences()[0][x].upper() == 'C':
                                state_set = ['C']
                            elif node.get_sequences()[0][x].upper() == 'G':
                                state_set = ['G']
                            elif node.get_sequences()[0][x].upper() == 'T':
                                state_set = ['T']
                            elif node.get_sequences()[0][x].upper() == 'R':
                                state_set = ['A','G']
                            elif node.get_sequences()[0][x].upper() == 'Y':
                                state_set = ['C','T']
                            elif node.get_sequences()[0][x].upper() == 'S':
                                state_set = ['C','G']
                            elif node.get_sequences()[0][x].upper() == 'W':
                                state_set = ['A','T']
                            elif node.get_sequences()[0][x].upper() == 'K':
                                state_set = ['G','T']
                            elif node.get_sequences()[0][x].upper() == 'M':
                                state_set = ['A','C']
                            elif node.get_sequences()[0][x].upper() == 'B':
                                state_set = ['C','G','T']
                            elif node.get_sequences()[0][x].upper() == 'D':
                                state_set = ['A','G','T']
                            elif node.get_sequences()[0][x].upper() == 'H':
                                state_set = ['A','C','T']
                            elif node.get_sequences()[0][x].upper() == 'V':
                                state_set = ['A','C','G']
                            elif node.get_sequences()[0][x].upper() == 'N':
                                state_set = ['A','C','G','T']
                            elif node.get_sequences()[0][x] == '?':
                                state_set = ['A','C','G','T']
                            elif node.get_sequences()[0][x] == '-':
                                state_set = ['A','C','G','T']
                            else:
                                print(node.get_sequences()[0][x])
                            node.set_state_set(x, state_set)
            invest_dist -= 1

        # Find one of the most parsimonious solutions (at random) to reconstruct
        # ancestral sequences (top down).
        # Deal with the root separately.
        for node in self.nodes:
            if node.get_is_root() == True:
                for x in range(seq_length):
                    state_set = node.get_state_set(x)
                    state = random.choice(state_set)
                    node.set_state(x, state)
                break
        # Deal with all other nodes (including terminal ones for consistency).
        invest_dist = 1
        while invest_dist <= self.max_dist_to_root:
            for node in self.nodes:
                if node.get_distance_to_root() == invest_dist:
                    for parent in self.nodes:
                        if node.get_parent_id() == parent.get_id():
                            for x in range(seq_length):
                                if parent.get_state(x) in node.get_state_set(x):
                                    node.set_state(x, parent.get_state(x))
                                else:
                                    node.set_state(x, random.choice(node.get_state_set(x)))
                            break
            invest_dist += 1
                    
        # Convert state sets to sequences for all nodes.
        for node in self.nodes:
            node.convert_states_to_sequence()

    def calculate_fitch_distances(self, transversions_only):
        for edge in self.edges:
            seqs = []
            for node in self.nodes:
                if node.get_id() in edge.get_node_ids():
                    seqs.append(node.get_sequences()[0])                    
            seq1 = XSeq(seqs[0])
            seq2 = XSeq(seqs[1])
            fitch_dist = seq1.get_distance_to(seq2, transversions_only)
            edge.set_fitch_distance(fitch_dist)

    def reduce(self, minimum_edge_length, minimum_node_size):
        # Unless an edge is removed already, check whether it should be removed
        # (this is the case if its Fitch distance is 0).
        for edge in self.edges:
            if edge.get_is_removed() == False:
                if edge.get_fitch_distance() < minimum_edge_length:
                    # Find the two nodes of that edge.
                    edge_nodes = []
                    for node in self.nodes:
                        if node.get_is_removed() == False:
                            if node.get_id() in edge.get_node_ids():
                                edge_nodes.append(node)
                    # Find out which of the two nodes of this edge is the parent
                    # of the other.
                    if edge_nodes[0].get_id() == edge_nodes[1].get_parent_id():
                        # edge_nodes[0] is always the parent of edge_nodes[1]
                        node_to_be_removed = edge_nodes[1]
                        node_to_be_kept = edge_nodes[0]
                    else:
                        # edge_nodes[1] is the parent of edge_nodes[0]
                        node_to_be_removed = edge_nodes[0]
                        node_to_be_kept = edge_nodes[1]
                    # Transfer properties from the node to be removed to the one to
                    # be kept.
                    node_to_be_kept.increase_size(node_to_be_removed.get_size())
                    for child_id in node_to_be_removed.get_child_ids():
                        node_to_be_kept.add_child_id(child_id)
                    for record_id in node_to_be_removed.get_record_ids():
                        node_to_be_kept.add_record_id(record_id)
                    for sequence in node_to_be_removed.get_sequences():
                        node_to_be_kept.add_sequence(sequence)
                    node_to_be_kept.remove_child_id(node_to_be_removed.get_id())
                    for pop in node_to_be_removed.get_pops():
                        node_to_be_kept.add_pop(pop)
                    # Find the children of node_to_be_removed and correct their parent_id.
                    for child in self.nodes:
                        if child.get_parent_id() == node_to_be_removed.get_id():
                            child.set_parent_id(node_to_be_removed.get_parent_id())
                    # Remove node_to_be_removed.
                    node_to_be_removed.set_is_removed(True)
                    # Find all edges that connected to the node, and replace their node id.
                    for downstream_edge in self.edges:
                        if downstream_edge != edge:
                            downstream_edge_node_ids = downstream_edge.get_node_ids()
                            if node_to_be_removed.get_id() in downstream_edge_node_ids:
                                downstream_edge_node_ids[0] = node_to_be_kept.get_id()
                                downstream_edge.set_node_ids(downstream_edge_node_ids)
                    # Remove this edge.
                    edge.set_is_removed(True)

        # Filter out the removed edges and nodes.
        new_edges = []
        for edge in self.edges:
            if edge.get_is_removed() == False:
                new_edges.append(edge)
        self.edges = new_edges
        new_nodes = []
        for node in self.nodes:
            if node.get_is_removed() == False:
                new_nodes.append(node)
        self.nodes = new_nodes

        # Recalculate the distances to the root.
        self.set_node_distances_to_root()

        # As a second reduction step, remove nodes (and their parental edges) if they
        # have no children, and a node size below the minimum node size. This is done
        # in a bottom-up way.
        if minimum_node_size > 1:
            invest_dist = self.max_dist_to_root
            while invest_dist >= 0:
                for node in self.nodes:
                    if node.get_distance_to_root() == invest_dist:
                        if node.get_child_ids() == [] and node.get_size() < minimum_node_size:
                            node.set_is_removed(True)
                            # Find the edge that connects to this node, and remove it.
                            for upstream_edge in self.edges:
                                upstream_edge_node_ids = upstream_edge.get_node_ids()
                                if node.get_id() in upstream_edge_node_ids:
                                    upstream_edge.set_is_removed(True)
                            # Find the parental node, and remove this child of it.
                            for parent in self.nodes:
                                if parent.get_id() == node.get_parent_id():
                                    parent.remove_child_id(node.get_id())
                invest_dist -= 1

            # Filter out the removed edges and nodes.
            new_edges = []
            for edge in self.edges:
                if edge.get_is_removed() == False:
                    new_edges.append(edge)
            self.edges = new_edges
            new_nodes = []
            for node in self.nodes:
                if node.get_is_removed() == False:
                    new_nodes.append(node)
            self.nodes = new_nodes

            # Recalculate the distances to the root.
            self.set_node_distances_to_root()

            # As a third reduction step, remove internal nodes with a size
            # below the minimum node size and exactly one downstream and one
            # upstream edge.
            change = True
            while change == True:
                change = False
                for node in self.nodes:
                    # All nodes except the root are considered.
                    if node.get_is_root() == False:
                        if node.get_size() < minimum_node_size and len(node.get_child_ids()) == 1:
                            node_id = node.get_id()
                            child_id = node.get_child_ids()[0]
                            parent_id = node.get_parent_id()
                            node.set_is_removed(True)
                            # Find the connecting edges.
                            for edge in self.edges:
                                if edge.get_is_removed() == False:
                                    edge_node_ids = edge.get_node_ids()
                                    if node_id in edge_node_ids:
                                        if child_id in edge_node_ids:
                                            downstream_edge_fitch_distance = edge.get_fitch_distance()
                                            # downstream_edge_length = edge.get_length()
                                            edge.set_is_removed(True)
                                        elif parent_id in edge_node_ids:
                                            upstream_edge_fitch_distance = edge.get_fitch_distance()
                                            # upstream_edge_length = edge.get_length()
                                            edge.set_is_removed(True)
                            new_edge = Edge("edge" + str(self.get_max_number_of_edges()) + "X")
                            self.increase_max_number_of_edges()
                            new_edge.set_node_ids([parent_id,child_id])
                            # new_edge_length = downstream_edge_length
                            # new_edge_length += upstream_edge_length
                            # new_edge.set_length(new_edge_length)
                            new_edge_fitch_distance = downstream_edge_fitch_distance 
                            new_edge_fitch_distance += upstream_edge_fitch_distance
                            new_edge.set_fitch_distance(new_edge_fitch_distance)
                            self.edges.append(new_edge)
                            # Find the parent and change its child ids.
                            for parent in self.nodes:
                                if parent.get_id() == parent_id:
                                    parent.remove_child_id(node_id)
                                    parent.add_child_id(child_id)
                                    break
                            # Find the child and replace its parent id.
                            for child in self.nodes:
                                if child.get_id() == child_id:
                                    child.set_parent_id(parent_id)
                            change = True
                            break

                if change == True:
                    # Filter out the removed edges and nodes.
                    new_edges = []
                    for edge in self.edges:
                        if edge.get_is_removed() == False:
                            new_edges.append(edge)
                    self.edges = new_edges
                    new_nodes = []
                    for node in self.nodes:
                        if node.get_is_removed() == False:
                            new_nodes.append(node)
                    self.nodes = new_nodes

            # In case the root node is below the minimum node size,
            # deal with it separately, as it was excluded above.
            change = False
            for old_root in self.nodes:
                if old_root.get_is_root() == True:
                    number_of_root_children = len(old_root.get_child_ids())
                    if old_root.get_size() < minimum_node_size and number_of_root_children == 2:
                        # Make one of the root children the new root, and make the other one
                        # the child of the new root.
                        old_root_child_ids = old_root.get_child_ids()
                        for old_root_child in self.nodes:
                            if old_root_child.get_id() == old_root_child_ids[0]:
                                new_root = old_root_child
                            elif old_root_child.get_id() == old_root_child_ids[1]:
                                new_root_child = old_root_child
                        old_root.set_is_root(False)
                        old_root.set_is_removed(True)
                        new_root.set_is_root(True)
                        new_root.set_parent_id('None')
                        new_root.add_child_id(new_root_child.get_id())
                        new_root_child.set_parent_id(new_root.get_id())
                        # Combine the two edges that connected to the old root to a new one.
                        # new_root_edge_length = 0
                        new_root_edge_fitch_distance = 0
                        for old_root_edge in self.edges:
                            old_root_edge_node_ids = old_root_edge.get_node_ids()
                            if old_root.get_id() in old_root_edge_node_ids:
                                # new_root_edge_length += old_root_edge.get_length()
                                new_root_edge_fitch_distance += old_root_edge.get_fitch_distance()
                                old_root_edge.set_is_removed(True)
                        new_root_edge = Edge("edge" + str(self.get_max_number_of_edges()) + "X")
                        self.increase_max_number_of_edges()
                        new_root_edge.set_node_ids([new_root.get_id(),new_root_child.get_id()])
                        # new_root_edge.set_length(new_root_edge_length)
                        new_root_edge.set_fitch_distance(new_root_edge_fitch_distance)
                        self.edges.append(new_root_edge)
                        change = True
                    break

            if change == True:
                # Filter out the removed edges and nodes.
                new_edges = []
                for edge in self.edges:
                    if edge.get_is_removed() == False:
                        new_edges.append(edge)
                self.edges = new_edges
                new_nodes = []
                for node in self.nodes:
                    if node.get_is_removed() == False:
                        new_nodes.append(node)
                self.nodes = new_nodes

            # Consider the case that the root is smaller than the minimum size,
            # but has only a single child.
            change = False
            for old_root in self.nodes:
                if old_root.get_is_root() == True:
                    number_of_root_children = len(old_root.get_child_ids())
                    if old_root.get_size() < minimum_node_size and number_of_root_children == 1:
                        # Make the only root child the new root, remove the root
                        # and the single root edge.
                        for new_root in self.nodes:
                            if new_root.get_id() == old_root.get_child_ids()[0]:
                                old_root.set_is_root(False)
                                old_root.set_is_removed(True)
                                new_root.set_is_root(True)
                                new_root.set_parent_id('None')
                        for old_root_edge in self.edges:
                            old_root_edge_node_ids = old_root_edge.get_node_ids()
                            if old_root.get_id() in old_root_edge_node_ids:
                                old_root_edge.set_is_removed(True)
                        change = True
                    break

            if change == True:
                # Filter out the removed edges and nodes.
                new_edges = []
                for edge in self.edges:
                    if edge.get_is_removed() == False:
                        new_edges.append(edge)
                self.edges = new_edges
                new_nodes = []
                for node in self.nodes:
                    if node.get_is_removed() == False:
                        new_nodes.append(node)
                self.nodes = new_nodes

        # Recalculate the distances to the root.
        self.set_node_distances_to_root()

        # Set the population sizes per population, for each node.
        for node in self.nodes:
            node.set_per_pop_sizes(self.pops)

    def position(self, algorithm, minimum_node_size, radius_multiplier):
        # Make a graph object from node and edge info.
        G = Graph()
        for node in self.nodes:
            G.add_node(node.get_id())
        for edge in self.edges:
            node_ids = edge.get_node_ids()
            G.add_edge(node_ids[0], node_ids[1])

        # Code below modified from networkx source.
        # This uses pygraphviz and neato to place nodes.
        A=pygraphviz.AGraph(name=G.name,strict=True,directed=False)
        # Default graph attributes
        A.graph_attr.update(G.graph.get('graph',{}))
        A.node_attr.update(G.graph.get('node',{}))
        A.edge_attr.update(G.graph.get('edge',{}))
        # Add nodes
        for n,nodedata in G.nodes(data=True):
            A.add_node(n,**nodedata)
        # Loop over edges
        for u,v,edgedata in G.edges_iter(data=True):
            str_edgedata=dict((k,str(v)) for k,v in edgedata.items())
            A.add_edge(u,v,**str_edgedata)
        A.layout(prog='neato')
        pos={}
        for n in G:
            node=pygraphviz.Node(A,n)
            try:
                xx,yy=node.attr["pos"].split(',')
                pos[n]=(float(xx),float(yy))
            except:
                print("no position for node",n)
                pos[n]=(0.0,0.0)
        
        # Define edges based on node positions.
        for edge in self.edges:
            edge_node_ids = edge.get_node_ids()
            for node_id, coordinates in pos.items():
                if node_id == edge_node_ids[0]:
                    x1 = coordinates[0]
                    y1 = coordinates[1]
                elif node_id == edge_node_ids[1]:
                    x2 = coordinates[0]
                    y2 = coordinates[1]
            uncorrected_delta_x = x2 - x1
            uncorrected_delta_y = y2 - y1
            inverse_scale_factor = math.sqrt(uncorrected_delta_x**2 + uncorrected_delta_y**2)
            unit_delta_x = uncorrected_delta_x/inverse_scale_factor
            unit_delta_y = uncorrected_delta_y/inverse_scale_factor
            edge.set_unit_delta_x(unit_delta_x)
            edge.set_unit_delta_y(unit_delta_y)

        # Get the sum of all fitch distances.
        sum_of_fitch_distances = 0
        for edge in self.edges:
            sum_of_fitch_distances += edge.get_fitch_distance()

        # Get the sum of all node sizes.
        sum_of_node_sizes = 0
        for node in self.nodes:
            sum_of_node_sizes += node.get_size()

        # Determine the scale factor of node size and edge length.
        if sum_of_fitch_distances == 0:
            scale_factor = 1
        else:
            scale_factor = sum_of_node_sizes/(2*sum_of_fitch_distances)

        # Set the position of the root node to 0,0.
        root_found = False
        for node in self.nodes:
            if node.get_is_root() == True:
                root_found = True
                node.set_x(0.0)
                node.set_y(0.0)
                break
        if root_found == False and len(self.nodes) > 0:
            print('WARNING: No root was found!')

        # Deal with all other nodes in a top-down sequence.
        invest_dist = 1
        while invest_dist <= self.max_dist_to_root:
            for node in self.nodes:
                if node.get_distance_to_root() == invest_dist:
                    node_id = node.get_id()
                    # Get the radius of this node.
                    node_radius = node.get_radius(minimum_node_size)*radius_multiplier
                    # Get the radius of the parent of this node:
                    parent_found = False
                    for parent in self.nodes:
                        if parent.get_id() == node.get_parent_id():
                            parent_found = True
                            parent_id = parent.get_id()
                            parent_radius = parent.get_radius(minimum_node_size)*radius_multiplier
                            parent_x = parent.get_x()
                            parent_y = parent.get_y()
                            break
                    if parent_found == False:
                        print('WARNING: Parent not found!')
                    # Get the edge that connects this node with its parent.
                    connecting_edge_found = False
                    for edge in self.edges:
                        edge_node_ids = edge.get_node_ids()
                        if edge_node_ids[0] == parent_id and edge_node_ids[1] == node_id:
                            connecting_edge_found = True
                            connecting_edge = edge
                            break
                    if connecting_edge_found == False:
                        print('WARNING: An edge was not found!')
                    connecting_edge_length = connecting_edge.get_fitch_distance() * scale_factor
                    total_edge_length = node_radius + parent_radius + connecting_edge_length
                    total_delta_x = total_edge_length * connecting_edge.get_unit_delta_x()
                    total_delta_y = total_edge_length * connecting_edge.get_unit_delta_y()
                    node.set_x(parent_x + total_delta_x)
                    node.set_y(parent_y + total_delta_y)
            invest_dist += 1
        self.is_positioned = True

    def to_svg(self, dim_x, dim_y, margin, minimum_node_size, radius_multiplier, colors, rest_color):
        if len(self.nodes) > 0:
            # Determine the two nodes with the greatest distance to each other.
            max_node_distance = 0
            index = 0
            if len(self.nodes) > 1:
                for node in self.nodes:
                    for other_node in self.nodes[(index+1):len(self.nodes)]:
                        delta_x = (node.get_x()-other_node.get_x())**2
                        delta_y = (node.get_y()-other_node.get_y())**2
                        node_distance = math.sqrt(delta_x + delta_y)
                        if node_distance > max_node_distance:
                            max_node_distance = node_distance
                            max_dist_node1 = node
                            max_dist_node2 = other_node
                if max_dist_node1.get_x() > max_dist_node2.get_x():
                    right_node = max_dist_node1
                    left_node = max_dist_node2
                else:
                    left_node = max_dist_node1
                    right_node = max_dist_node2
                # Determine the maximum horizontal and vertical distance between these two nodes.
                max_distance_delta_x = right_node.get_x() - left_node.get_x()
                max_distance_delta_y = right_node.get_y() - left_node.get_y()
                # Find the center between these two nodes, it becomes the center of rotation.
                rotation_center_x = 0.5 * (left_node.get_x() + right_node.get_x())
                rotation_center_y = 0.5 * (left_node.get_y() + right_node.get_y())
                # The angle of the following rotation is determined.
                # The following two cases are treated seperately in order to avoid division by zero.
                # Find the current (unscaled) min and max of x and y.
                # In the last case, the arctan helps to find the rotation angle.
                if max_distance_delta_y == 0.0:
                    w = 0
                elif max_distance_delta_x == 0.0:
                    w = 0.5 * math.pi
                elif max_distance_delta_y < 0.0:
                    w = math.atan(-max_distance_delta_y/max_distance_delta_x)
                elif max_distance_delta_y > 0.0:
                    w = - math.atan(max_distance_delta_y/max_distance_delta_x)
                # Finally, the rotation is performed on all nodes.
                for node in self.nodes:
                    node.rotate(w, rotation_center_x, rotation_center_y)
            # After the rotation, get once more the dimensions of the graph.
            unscaled_min_x = 0
            unscaled_max_x = 0
            unscaled_min_y = 0
            unscaled_max_y = 0
            for node in self.nodes:
                if node.get_x() - node.get_radius(minimum_node_size)*radius_multiplier < unscaled_min_x:
                    unscaled_min_x = node.get_x() - node.get_radius(minimum_node_size)*radius_multiplier
                if node.get_x() + node.get_radius(minimum_node_size)*radius_multiplier > unscaled_max_x:
                    unscaled_max_x = node.get_x() + node.get_radius(minimum_node_size)*radius_multiplier
                if node.get_y() - node.get_radius(minimum_node_size)*radius_multiplier < unscaled_min_y:
                    unscaled_min_y = node.get_y() - node.get_radius(minimum_node_size)*radius_multiplier
                if node.get_y() + node.get_radius(minimum_node_size)*radius_multiplier > unscaled_max_y:
                    unscaled_max_y = node.get_y() + node.get_radius(minimum_node_size)*radius_multiplier
            # Determine the scale factor for all x and y values.
            scale_factor_x = (dim_x - 2*margin)/(unscaled_max_x - unscaled_min_x)
            scale_factor_y = (dim_y - 2*margin)/(unscaled_max_y - unscaled_min_y)
            scale_factor = scale_factor_x
            if scale_factor_y < scale_factor:
                scale_factor = scale_factor_y
            scaled_range_x = scale_factor*(unscaled_max_x - unscaled_min_x)
            scaled_range_y = scale_factor*(unscaled_max_y - unscaled_min_y)
            # Determine the extra margin on the left or on the top (one of them is 0).
            extra_x_margin = ((dim_x - 2*margin) - scaled_range_x)/2
            extra_y_margin = ((dim_y - 2*margin) - scaled_range_y)/2
            # Reposition all nodes based on extra margin and scaling factor.
            for node in self.nodes:
                node_x = node.get_x() - unscaled_min_x
                node_y = node.get_y() - unscaled_min_y
                node_radius = node.get_radius(minimum_node_size)
                node.set_x(margin + extra_x_margin + scale_factor * node_x)
                node.set_y(margin + extra_y_margin + scale_factor * node_y)
                node.set_radius(node_radius * scale_factor)
            # Check once more the top and bottom margins, now after all scaling:
            top_margin = dim_y
            bottom_margin = dim_y
            for node in self.nodes:
                if node.get_y() - node.get_radius(minimum_node_size)*radius_multiplier < top_margin:
                    top_margin = node.get_y() - node.get_radius(minimum_node_size)*radius_multiplier
                if dim_y - (node.get_y() + node.get_radius(minimum_node_size)*radius_multiplier) < bottom_margin:
                    bottom_margin = dim_y - (node.get_y() + node.get_radius(minimum_node_size)*radius_multiplier)
            adjusted_top_margin = top_margin % 10 + 10
            adjusted_bottom_margin = bottom_margin % 10 + 10
            top_margin_to_be_removed = top_margin - adjusted_top_margin
            bottom_margin_to_be_removed = bottom_margin - adjusted_bottom_margin
            dim_y = int(dim_y - (top_margin_to_be_removed + bottom_margin_to_be_removed))
            for node in self.nodes:
                node_y = node.get_y()
                node.set_y(node_y - top_margin_to_be_removed)

            # Definitions for the svg.
            stroke_width = 1.0
            stroke_color = '93a1a1' # base1
            font_color = '002b36'  # base03
            font_unit = 10
            font = '\'Helvetica\''
            text_y_correction = 0.35
            dot_radius = 1.5

            # Start writing the svg string.
            svg_string = ''
            svg_string += '<!-- The haplotype genealogy graph: info start -->\n'
            svg_string += '<!-- Populations:\n'
            for x in range(0,len(pops)):
                if x < len(colors):
                    svg_string += pops[x] + ': #' + colors[x] + '\n'
                else:
                    svg_string += pops[x] + ': #' + rest_color + '\n'
            svg_string += '-->\n'
            svg_string += '<!-- The haplotype genealogy graph: info end -->\n'
            svg_string += '<!-- The haplotype genealogy graph: SVG string start -->\n'
            svg_string += '<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"'
            svg_string += str(dim_x)
            svg_string += '\" height=\"'
            svg_string += str(dim_y)
            svg_string += '\" viewBox=\"'
            svg_string += '0 0 '
            svg_string += str(dim_x)
            svg_string += ' '
            svg_string += str(dim_y)
            svg_string += '\" xmlns:xlink=\"htttp://www.w3.org/1999/xlink\">\n\n'
            svg_string += '  <defs>\n'
            svg_string += '    <style type=\"text/css\">\n'
            svg_string += '      <![CDATA[\n'
            svg_string += '        text {font-weight:normal; fill:#'
            svg_string += font_color
            svg_string += ';}\n'
            svg_string += '        line {fill:none; stroke:#'
            svg_string += stroke_color
            svg_string += '; stroke-width:'
            svg_string += str(stroke_width)
            svg_string += 'px;}\n'
            svg_string += '        circle {stroke:#'
            svg_string += stroke_color
            svg_string += '; stroke-width:'
            svg_string += str(stroke_width)
            svg_string += 'px;}\n'
            svg_string += '      ]]>\n'
            svg_string += '    </style>\n\n'
            svg_string += '    <radialGradient id=\"radgrad\"\n'
            svg_string += '      cx=\".9\" cy=\".1\" r=\"1\">\n'
            svg_string += '      <stop  offset=\"0\" style=\"stop-color:grey\"/>\n'
            svg_string += '      <stop  offset=\"1\" style=\"stop-color:black\"/>\n'
            svg_string += '    </radialGradient>\n\n'
            svg_string += '    <!--Masks-->\n'
            count = 0
            for node in self.nodes:
                svg_string += '    <mask id=\"m_'
                svg_string += str(count)
                svg_string += '\"><circle fill=\"url(#radgrad)\" cx=\"'
                svg_string += str(node.get_x())
                svg_string += '\" cy=\"'
                svg_string += str(node.get_y())
                svg_string += '\" r=\"'
                svg_string += str(node.get_radius(minimum_node_size)*radius_multiplier)
                svg_string += '\"/></mask>\n'
                count += 1
            svg_string += '  </defs>\n\n'

            # Add edges.
            svg_string += '  <!--Edges-->\n'
            for edge in self.edges:
                edge_node_ids = edge.get_node_ids()
                for node in self.nodes:
                    if node.get_id() == edge_node_ids[0]:
                        x1 = node.get_x()
                        y1 = node.get_y()
                    elif node.get_id() == edge_node_ids[1]:
                        x2 = node.get_x()
                        y2 = node.get_y()
                svg_string += '  <line x1=\"'
                svg_string += str(x1)
                svg_string += '\" y1=\"'
                svg_string += str(y1)
                svg_string += '\" x2=\"'
                svg_string += str(x2)
                svg_string += '\" y2=\"'
                svg_string += str(y2)
                svg_string += '\"/>\n'
            svg_string += '\n'

            # Add dots to edges to mark substitutions.
            svg_string += '  <!--Dots-->\n'
            for edge in self.edges:
                if edge.get_fitch_distance() > 1:
                    edge_node_ids = edge.get_node_ids()
                    for node in self.nodes:
                        if node.get_id() == edge_node_ids[0]:
                            node1_x = node.get_x()
                            node1_y = node.get_y()
                            node1_r = node.get_radius(minimum_node_size)*radius_multiplier
                        elif node.get_id() == edge_node_ids[1]:
                            node2_x = node.get_x()
                            node2_y = node.get_y()
                            node2_r = node.get_radius(minimum_node_size)*radius_multiplier
                    if node2_x > node1_x:
                        left_node_x = node1_x
                        left_node_y = node1_y
                        left_node_r = node1_r
                        right_node_x = node2_x
                        right_node_y = node2_y
                        right_node_r = node2_r
                    else:
                        right_node_x = node1_x
                        right_node_y = node1_y
                        right_node_r = node1_r
                        left_node_x = node2_x
                        left_node_y = node2_y
                        left_node_r = node2_r
                    edge_delta_x = right_node_x - left_node_x
                    edge_delta_y = right_node_y - left_node_y
                    edge_total_length = math.sqrt(edge_delta_x**2 + edge_delta_y**2)
                    edge_unit_delta_x = edge_delta_x / edge_total_length
                    edge_unit_delta_y = edge_delta_y / edge_total_length
                    if edge_delta_y == 0: # edge is horizontal.
                        start_x = left_node_x + left_node_r
                        start_y = left_node_y
                        end_x = right_node_x - right_node_r
                        end_y = right_node_y
                    elif edge_delta_y > 0: # edge goes down.
                        if edge_delta_x == 0: # edge goes down vertically.
                            start_x = left_node_x
                            start_y = left_node_y + left_node_r
                            end_x = right_node_x
                            end_y = right_node_x - right_node_r
                        else: # edge goes down, but not vertically.
                            start_x = left_node_x + edge_unit_delta_x * left_node_r
                            start_y = left_node_y + edge_unit_delta_y * left_node_r
                            end_x = right_node_x - edge_unit_delta_x * right_node_r
                            end_y = right_node_y - edge_unit_delta_y * right_node_r
                    else: # edge goes up.
                        if edge_delta_x == 0: # edge goes up vertically.
                            start_x = left_node_x
                            start_y = left_node_y - left_node_r
                            end_x = right_node_x
                            end_y = right_node_x + right_node_r
                        else: # edge goes up, but not vertically.
                            start_x = left_node_x + edge_unit_delta_x * left_node_r
                            start_y = left_node_y + edge_unit_delta_y * left_node_r
                            end_x = right_node_x - edge_unit_delta_x * right_node_r
                            end_y = right_node_y - edge_unit_delta_y * right_node_r
                    substitution_x = (end_x - start_x)/edge.get_fitch_distance()
                    substitution_y = (end_y - start_y)/edge.get_fitch_distance()
                    for i in range(edge.get_fitch_distance()-1):
                        dot_x = start_x + (i+1) * substitution_x
                        dot_y = start_y + (i+1) * substitution_y
                        svg_string += '  <circle fill=\"#'
                        svg_string += stroke_color
                        svg_string += '\" cx=\"'
                        svg_string += str(dot_x)
                        svg_string += '\" cy=\"'
                        svg_string += str(dot_y)
                        svg_string += '\" r=\"'
                        svg_string += str(dot_radius)
                        svg_string += '\"/>\n'
            svg_string += '\n'

            # Add the nodes.
            svg_string += '  <!--Nodes-->\n'
            count = 0
            for node in self.nodes:
                if node.get_size() > 0:
                    svg_string += '  <!--Haplotype ID: '
                    svg_string += str(count+1)
                    svg_string += '-->\n'
                    svg_string += '  <circle fill=\"#'
                    svg_string += rest_color
                    svg_string += '\" cx=\"'
                    svg_string += str(node.get_x())
                    svg_string += '\" cy=\"'
                    svg_string += str(node.get_y())
                    svg_string += '\" r=\"'
                    svg_string += str(node.get_radius(minimum_node_size)*radius_multiplier)
                    svg_string += '\"/>\n'
                    old_p =  0
                    new_p = 0
                    for s in range(len(self.pops)):
                        proportion = node.get_per_pop_sizes()[s] / node.get_size()
                        if proportion > 0 and proportion < 1:
                            new_p = old_p + proportion
                            pie_start_x = node.get_x() + (node.get_radius(minimum_node_size)*radius_multiplier - 0.5*stroke_width)*math.sin(old_p*2*math.pi)
                            pie_start_y = node.get_y() - (node.get_radius(minimum_node_size)*radius_multiplier - 0.5*stroke_width)*math.cos(old_p*2*math.pi)
                            pie_stop_x = node.get_x() + (node.get_radius(minimum_node_size)*radius_multiplier - 0.5*stroke_width)*math.sin(new_p*2*math.pi)
                            pie_stop_y = node.get_y() - (node.get_radius(minimum_node_size)*radius_multiplier - 0.5*stroke_width)*math.cos(new_p*2*math.pi)
                            svg_string += '  <path fill=\"#'
                            if s >= len(colors):
                                svg_string += rest_color
                            else:
                                svg_string += colors[s]
                            svg_string += '\" d=\"M '
                            svg_string += str(node.get_x())
                            svg_string += ','
                            svg_string += str(node.get_y())
                            svg_string += ' L '
                            svg_string += str(pie_start_x)
                            svg_string += ','
                            svg_string += str(pie_start_y)
                            svg_string += ' A '
                            svg_string += str(node.get_radius(minimum_node_size)*radius_multiplier - 0.5*stroke_width)
                            svg_string += ' '
                            svg_string += str(node.get_radius(minimum_node_size)*radius_multiplier - 0.5*stroke_width)
                            svg_string += ' 0 '
                            svg_string += str(round(proportion))
                            svg_string += ' 1 '
                            svg_string += str(pie_stop_x)
                            svg_string += ', '
                            svg_string += str(pie_stop_y)
                            svg_string += ' z\"/>\n'
                            old_p = new_p
                        elif proportion == 1.0:
                            svg_string += '  <circle fill=\"#'
                            if s > len(colors)-1:
                                svg_string += rest_color
                            else:
                                svg_string += colors[s]
                            svg_string += '\" cx=\"'
                            svg_string += str(node.get_x())
                            svg_string += '\" cy=\"'
                            svg_string += str(node.get_y())
                            svg_string += '\" r=\"'
                            svg_string += str(node.get_radius(minimum_node_size)*radius_multiplier)
                            svg_string += '\"/>\n'
                        elif proportion > 1 or proportion < 0:
                            warn_string = ''
                            warn_string += 'WARNING: A proportion was found to '
                            warn_string += 'be less than 0 or greater than 1!'
                            print(warn_string)
                    svg_string += '\n'
                else:
                    svg_string += '  <circle fill=\"#'
                    svg_string += stroke_color
                    svg_string += '\" cx=\"'
                    svg_string += str(node.get_x())
                    svg_string += '\" cy=\"'
                    svg_string += str(node.get_y())
                    svg_string += '\" r=\"'
                    svg_string += str(dot_radius)
                    svg_string += '\"/>\n'
                count += 1

            # Add the gradients.
            svg_string += '  <!--Gradients-->\n'
            count = 0
            for node in self.nodes:
                svg_string += '  <circle fill=\"white\" cx=\"'
                svg_string += str(node.get_x())
                svg_string += '\" cy=\"'
                svg_string += str(node.get_y())
                svg_string += '\" r=\"'
                svg_string += str(node.get_radius(minimum_node_size)*radius_multiplier)
                svg_string += '\" mask=\"url(#m_'
                svg_string += str(count)
                svg_string += ')\"/>\n'
                count += 1
            svg_string += '\n'

            # Add the node labels.
            svg_string += '  <!--Node labels-->\n'
            count = 0
            for node in self.nodes:
                if count < 9:
                    text_x_correction = -0.2
                elif count == 10:
                    text_x_correction = -0.52
                else:
                    text_x_correction = -0.60
                svg_string += '  <text x=\"'
                svg_string += str(node.get_x() + font_unit * text_x_correction)
                svg_string += '\" y=\"'
                svg_string += str(node.get_y() + font_unit * text_y_correction)
                svg_string += '\" font-family=\"'
                svg_string += str(font)
                svg_string += '\" font-size=\"'
                svg_string += str(font_unit)
                svg_string += 'px\">'
                svg_string += str(count+1)
                svg_string += '</text>\n'
                count += 1
            svg_string += '\n'

            # Finish and return the svg string.
            svg_string += '</svg>\n'
            svg_string += '<!-- The haplotype genealogy graph: SVG string end -->\n'
            return svg_string

        # If not a single node exists (cause they all fall below the size limit).
        else:
            return 'No nodes found with sufficient size for display!'

    def assign_progeny_ids(self):
        for node in self.nodes:
            node_progeny_ids = []
            node_ids_to_follow = []
            for child_id in node.get_child_ids():
                node_progeny_ids.append(child_id)
                node_ids_to_follow.append(child_id)
            while len(node_ids_to_follow) > 0:
                for node_id_to_follow in node_ids_to_follow:
                    node_ids_to_follow.remove(node_id_to_follow)
                    for progeny_node in self.nodes:
                        if progeny_node.get_id() == node_id_to_follow:
                            for child_id in progeny_node.get_child_ids():
                                node_progeny_ids.append(child_id)
                                node_ids_to_follow.append(child_id)
            node.set_progeny_ids(node_progeny_ids)

    def get_number_of_internal_nodes(self):
        number_of_internal_nodes = 0
        for node in self.nodes:
            if 'internalNode' in node.get_id():
                number_of_internal_nodes += 1
        return number_of_internal_nodes

    def get_gsi(self, pop):
        pop_terminals = []
        for node in self.nodes:
            if pop in node.get_id():
                pop_terminals.append(node.get_id())
        gsis = []
        if len(pop_terminals) == 0:
            return None
        if len(pop_terminals) == 1:
            return 1
        else:
            for node in self.nodes:
                # print(node.info())
                if 'internalNode' in node.get_id():
                    if node.extant_progeny_ids_contain_all_of(pop_terminals):
                        number_of_nodes_actually_required = 1
                        for progeny_id in node.get_progeny_ids():
                            for progeny_node in self.nodes:
                                if progeny_node.get_id() == progeny_id:
                                    if progeny_node.extant_progeny_ids_contain_any_of(pop_terminals):
                                        number_of_nodes_actually_required += 1
                    elif node.extant_progeny_ids_contain_none_of(pop_terminals):
                        number_of_nodes_actually_required = 0
                        # Find internal nodes that are above this node (but exclude the root).
                        for anti_progeny_node in self.nodes:
                            if 'internalNode' in anti_progeny_node.get_id():
                                if anti_progeny_node.get_id() not in node.get_progeny_ids():
                                    if anti_progeny_node.get_is_root() == False:
                                        if anti_progeny_node.extant_progeny_ids_contain_any_of(pop_terminals):
                                            number_of_nodes_actually_required += 1
                        # Check whether the root is also required.
                        number_of_root_children_with_any_pop_terminals = 0
                        for root in self.nodes:
                            if root.get_is_root():
                                for root_child_id in root.get_child_ids():
                                    for child in self.nodes:
                                        if child.get_id() == root_child_id:
                                            if 'internalNode' in child.get_id():
                                                if child.extant_progeny_ids_contain_any_of(pop_terminals):
                                                    number_of_root_children_with_any_pop_terminals += 1
                                            else:
                                                if pop in child.get_id():
                                                    number_of_root_children_with_any_pop_terminals += 1
                                            break
                                break
                        if number_of_root_children_with_any_pop_terminals > 1:
                            number_of_nodes_actually_required += 1
                    else:
                        number_of_nodes_actually_required = None
                    if number_of_nodes_actually_required != None:
                        number_of_nodes_minimally_required = len(pop_terminals) - 1
                        # The + 1 in the next line is needed since the root is the only node that
                        # has three children (it's a polytomy).
                        number_of_nodes_maximally_required = self.get_number_of_internal_nodes() + 1
                        obs_gs = number_of_nodes_minimally_required/number_of_nodes_actually_required
                        max_gs = 1
                        min_gs = number_of_nodes_minimally_required/number_of_nodes_maximally_required
                        gsi = (obs_gs - min_gs) / (max_gs - min_gs)
                        gsis.append(gsi)
            max_gsi = None
            for gsi in gsis:
                if gsi != None:
                    if max_gsi == None:
                        max_gsi = 0
                    if gsi > max_gsi:
                        max_gsi = gsi
            return max_gsi

    def info(self):
        info_string = ''
        info_string += 'Tree'.ljust(20)
        info_string += '\n'
        info_string += 'Number of nodes:'.ljust(20)
        info_string += str(self.get_number_of_nodes())
        info_string += '\n'
        info_string += 'Number of edges:'.ljust(20)
        info_string += str(self.get_number_of_edges())
        info_string += '\n'
        return info_string

    def check_edges(self):
        all_sorted = True
        for edge in self.edges:
            edge_node_ids = edge.get_node_ids()
            for node in self.nodes:
                if node.get_id() == edge_node_ids[0]:
                    dist1 = node.get_distance_to_root()
                elif node.get_id() == edge_node_ids[1]:
                    dist2 = node.get_distance_to_root()
            if dist1 > dist2:
                all_sorted = False
        return all_sorted


# The Node class.
class Node(object):

    def __init__(self, id, is_root, pops):
        self.id = id
        self.is_root = is_root
        self.sequences = []
        self.child_ids = []
        self.number_of_children = 0
        self.state_sets = []
        self.states = []
        self.size = 0
        self.record_ids = []
        self.radius = 0
        self.is_removed = False
        self.pops = []
        self.per_pop_sizes = []
        self.progeny_ids = []
        self.distance_to_root = None
        for pop in pops:
            if pop in self.id:
                self.pops.append(pop)
        if self.pops == [] and self.id[:12] != 'internalNode':
            self.pops.append('unknown')
        self.x = 'None'
        self.y = 'None'

    def get_id(self):
        return self.id

    def set_size(self, size):
        self.size = size
        self.set_radius(0.5*math.sqrt(self.size))

    def get_size(self):
        return self.size

    def increase_size(self, increase_size):
        self.size += increase_size
        self.set_radius(0.5*math.sqrt(self.size))

    def add_record_id(self, record_id):
        self.record_ids.append(record_id)

    def get_record_ids(self):
        return self.record_ids

    def set_radius(self, radius):
        self.radius = radius

    def get_radius(self, minimum_node_size):
        if self.size >= minimum_node_size:
            return self.radius
        else:
            return 0

    def set_x(self, x):
        self.x = x

    def get_x(self):
        return self.x

    def set_y(self, y):
        self.y = y

    def get_y(self):
        return self.y

    def switch_x_and_y(self):
        self.x, self.y = self.y, self.x

    def rotate(self, w, rotation_center_x, rotation_center_y):
        squared_distance_to_center_x = (self.x - rotation_center_x)**2
        squared_distance_to_center_y = (self.y - rotation_center_y)**2
        distance_to_center = math.sqrt(squared_distance_to_center_x + squared_distance_to_center_y)
        old_delta_x = self.x - rotation_center_x
        old_delta_y = self.y - rotation_center_y
        if old_delta_y >= 0:
            new_delta_x = distance_to_center * (math.sin(math.asin(old_delta_x/distance_to_center)-w))
            self.x = new_delta_x + rotation_center_x
        else:
            new_delta_x = distance_to_center * (-math.sin(math.asin(-old_delta_x/distance_to_center)-w))
            self.x = new_delta_x + rotation_center_x
        if old_delta_x <= 0:
            new_delta_y = distance_to_center * (math.sin(math.asin(old_delta_y/distance_to_center)-w))
            self.y = new_delta_y + rotation_center_y
        else:
            new_delta_y = distance_to_center * (-math.sin(math.asin(-old_delta_y/distance_to_center)-w))
            self.y = new_delta_y + rotation_center_y

    def set_sequences(self, sequences):
        self.sequences = sequences

    def add_sequence(self, sequence):
        self.sequences.append(sequence)

    def get_sequences(self):
        return self.sequences

    def set_parent_id(self, parent_id):
        self.parent_id = parent_id

    def get_parent_id(self):
        return self.parent_id

    def add_child_id(self, child_id):
        self.number_of_children += 1
        self.child_ids.append(child_id)

    def remove_child_id(self, child_id):
        self.number_of_children -= 1
        self.child_ids.remove(child_id)

    def get_child_ids(self):
        return self.child_ids

    def get_number_of_children(self):
        return self.number_of_children

    def set_pops(self, pops):
        self.pops = pops

    def get_pops(self):
        return self.pops

    def add_pop(self, pop):
        self.pops.append(pop)

    def set_per_pop_sizes(self, pops):
        pop_pos = 0
        for pop in pops:
            this_pop_count = 0
            for own_pop in self.pops:
                if own_pop == pop:
                    this_pop_count += 1
            self.per_pop_sizes.append(this_pop_count)
            pop_pos += 1

    def get_per_pop_sizes(self):
        return self.per_pop_sizes

    def set_distance_to_root(self, distance_to_root):
        self.distance_to_root = distance_to_root

    def get_distance_to_root(self):
        return self.distance_to_root

    def set_is_root(self, is_root):
        self.is_root = is_root

    def get_is_root(self):
        if self.is_root == True:
            return True
        else:
            return False

    def prepare_state_sets(self, seq_length):
        for x in range(seq_length):
            self.state_sets.append(['None'])

    def set_state_set(self, pos, state):
        self.state_sets[pos] = state

    def get_state_set(self, pos):
        return self.state_sets[pos]

    def prepare_states(self, seq_length):
        for x in range(seq_length):
            self.states.append('None')

    def set_state(self, pos, state):
        self.states[pos] = state

    def get_state(self, pos):
        return self.states[pos]

    def convert_states_to_sequence(self):
        seq_string = ''
        for state in self.states:
            seq_string += state
        if self.sequences == []:
            self.sequences.append(seq_string)
        else:
            self.sequences[0] = seq_string

    def set_is_removed(self, is_removed):
        if is_removed == True:
            self.is_removed = True
        elif is_removed == False:
            self.is_removed = False
        else:
            print('WARNING: Unexpected value for boolean variable.')

    def get_is_removed(self):
        return self.is_removed

    def set_progeny_ids(self, progeny_ids):
        self.progeny_ids = progeny_ids

    def get_progeny_ids(self):
        return self.progeny_ids

    def get_extant_progeny_ids(self):
        if self.progeny_ids == []:
            return []
        else:
            extant_progeny_ids = []
            for progeny_id in self.progeny_ids:
                if 'internalNode' not in progeny_id:
                    extant_progeny_ids.append(progeny_id)
            return extant_progeny_ids

    def extant_progeny_ids_contain_all_of(self, terminal_ids):
        extant_progeny_ids_contain_all = True
        for terminal_id in terminal_ids:
            if terminal_id not in self.get_extant_progeny_ids():
                extant_progeny_ids_contain_all = False
        return extant_progeny_ids_contain_all

    def extant_progeny_ids_contain_any_of(self, terminal_ids):
        extant_progeny_ids_contain_any = False
        for terminal_id in terminal_ids:
            if terminal_id in self.get_extant_progeny_ids():
                extant_progeny_ids_contain_any = True
        return extant_progeny_ids_contain_any

    def extant_progeny_ids_contain_none_of(self, terminal_ids):
        extant_progeny_ids_contain_none = True
        for terminal_id in terminal_ids:
            if terminal_id in self.get_extant_progeny_ids():
                extant_progeny_ids_contain_none = False
        return extant_progeny_ids_contain_none

    def info(self):
        info_string = ''
        info_string += 'Node id:'.ljust(20)
        info_string += self.id
        info_string += '\n'
        info_string += 'Root:'.ljust(20)
        info_string += str(self.is_root)
        info_string += '\n'
        info_string += 'Parent id:'.ljust(20)
        info_string += str(self.parent_id)
        info_string += '\n'
        info_string += 'Child ids:'.ljust(20)
        if len(self.child_ids) > 0:
            for child_id in self.child_ids:
                info_string += child_id
                info_string += ', '
            info_string = info_string[:-2]
        else:
            info_string += 'None'
        info_string += '\n'
        info_string += 'Extant progeny ids:'.ljust(20)
        extant_progeny_ids = self.get_extant_progeny_ids()
        if len(extant_progeny_ids) > 0:
            for extant_progeny_id in extant_progeny_ids:
                info_string += extant_progeny_id
                info_string += ', '
            info_string = info_string[:-2]
        else:
            info_string += 'None'
        info_string += '\n'
        info_string += 'Distance to root:'.ljust(20)
        info_string += str(self.distance_to_root)
        info_string += '\n'
        info_string += 'Size:'.ljust(20)
        info_string += str(self.size)
        info_string += '\n'
        info_string += 'Radius:'.ljust(20)
        info_string += str(self.radius)
        info_string += '\n'
        info_string += 'Pops:'.ljust(20)
        if len(self.pops) > 0:
            for pop in self.pops:
                info_string += pop
                info_string += ', '
            info_string = info_string[:-2]
        else:
            info_string += 'None'
        info_string += '\n'
        info_string += 'Pop sizes:'.ljust(20)
        if len(self.per_pop_sizes) > 0:
            for pop_size in self.per_pop_sizes:
                info_string += str(pop_size)
                info_string += ', '
            info_string = info_string[:-2]
        else:
            info_string += 'None'
        info_string += '\n'
        info_string += 'Sequences:'.ljust(20)
        if self.sequences != []:
            sequence_string = ''
            for sequence in self.sequences:
                sequence_string += sequence + ', '
            sequence_string = sequence_string[:-2]
            info_string += sequence_string
        else:
            info_string += 'None'
        info_string += '\n'
        # info_string += 'State sets:'.ljust(20)
        # for state_set in self.state_sets:
        #     info_string += '['
        #     for state in state_set:
        #         info_string += state
        #         info_string += ','
        #     info_string = info_string[:-1]
        #     info_string += ']'
        # info_string += '\n'
        info_string += 'Coordinates:'.ljust(20)
        info_string += str(self.x)
        info_string += ', '
        info_string += str(self.y)
        info_string += '\n'
        return info_string


# The Edge class.
class Edge(object):

    def __init__(self, id):
        self.id = id
        self.node_ids = []
        self.is_removed = False
        self.unit_delta_x = None
        self.unit_delta_y = None
        self.fitch_distance = None

    def get_id(self):
        return self.id

    def set_node_ids(self, node_ids):
        self.node_ids = node_ids

    def get_node_ids(self):
        return self.node_ids

    # def set_length(self, length):
    #     self.length = length

    # def get_length(self):
    #     return self.length

    def set_fitch_distance(self, fitch_dist):
        self.fitch_distance = fitch_dist

    def get_fitch_distance(self):
        return self.fitch_distance

    def set_unit_delta_x(self, unit_delta_x):
        self.unit_delta_x = unit_delta_x

    def get_unit_delta_x(self):
        return self.unit_delta_x

    def set_unit_delta_y(self, unit_delta_y):
        self.unit_delta_y = unit_delta_y

    def get_unit_delta_y(self):
        return self.unit_delta_y

    def set_is_removed(self, is_removed):
        if is_removed == True:
            self.is_removed = True
        elif is_removed == False:
            self.is_removed = False
        else:
            print('WARNING: Unexpected value for boolean variable.')

    def get_is_removed(self):
        return self.is_removed

    def info(self):
        info_string = ''
        info_string += 'Edge id:'.ljust(20)
        info_string += self.id
        info_string += '\n'
        info_string += 'Edge node 1 id:'.ljust(20)
        info_string += self.node_ids[0]
        info_string += '\n'
        info_string += 'Edge node 2 id:'.ljust(20)
        info_string += self.node_ids[1]
        info_string += '\n'
        info_string += 'Edge length:'.ljust(20)
        info_string += str(self.length)
        info_string += '\n'
        info_string += 'Unit delta x:'.ljust(20)
        info_string += str(self.unit_delta_x)
        info_string += '\n'
        info_string += 'Unit delta y:'.ljust(20)
        info_string += str(self.unit_delta_y)
        info_string += '\n'
        if self.fitch_distance != None:
            info_string += 'Fitch distance:'.ljust(20)
            info_string += str(self.fitch_distance)
            info_string += '\n'
        return info_string


# Expand the class Seq to include a distance measure method.
class XSeq(Seq):

    def get_distance_to(self, seq, transversions_only):
        distance = 0
        if len(self) != len(seq):
            print("ERROR: Both sequences must be of the same length!")
            sys.exit(1)
        for x in range(0,len(self)):
            if transversions_only:
                if self[x] in ['A', 'G', 'R']:
                    if seq[x] in ['A', 'G', 'S', 'K', 'B', 'W', 'R', 'M', 'D', 'H', 'V', 'N', '-', '?']:
                        distance += 0
                    elif seq[x] in ['C', 'T', 'Y']:
                        distance += 1
                    else:
                        print("Found unexpected base!")
                elif self[x] in ['C', 'T', 'Y']:
                    if seq[x] in ['C', 'T', 'Y', 'S', 'K', 'B', 'W', 'M', 'D', 'H', 'V', 'N', '-', '?']:
                        distance += 0
                    elif seq[x] in ['A', 'G', 'R']:
                        distance += 1
                    else:
                        print("Found unexpected base!")
                elif self[x] in ['S', 'K', 'B', 'W', 'M', 'D', 'H', 'V', 'N', '-', '?']:
                    distance += 0
                else:
                    print("Found unexpected base!")
            else:
                if self[x] == 'A':
                    if seq[x] == 'A':
                        distance += 0
                    elif seq[x] in ['C', 'G', 'T', 'Y', 'S', 'K', 'B']:
                        distance += 1
                    elif seq[x] in ['W', 'R', 'M']:
                        distance += 1/2
                    elif seq[x] in ['D', 'H', 'V']:
                        distance += 2/3
                    elif seq[x] in ['N', '-', '?']:
                        distance += 3/4
                    else:
                        print("Found unexpected base!")
                elif self[x] == 'C':
                    if seq[x] == 'C':
                        distance += 0
                    elif seq[x] in ['A', 'G', 'T', 'R', 'W', 'K', 'D']:
                        distance += 1
                    elif seq[x] in ['Y', 'S', 'M']:
                        distance += 1/2
                    elif seq[x] in ['B', 'H', 'V']:
                        distance += 2/3
                    elif seq[x] in ['N', '-', '?']:
                        distance += 3/4
                    else:
                        print("Found unexpected base!")
                elif self[x] == 'G':
                    if seq[x] == 'G':
                        distance += 0
                    elif seq[x] in ['A', 'C', 'T', 'Y', 'W', 'M', 'H']:
                        distance += 1
                    elif seq[x] in ['R', 'S', 'K']:
                        distance += 1/2
                    elif seq[x] in ['B', 'D', 'V']:
                        distance += 2/3
                    elif seq[x] in ['N', '-', '?']:
                        distance += 3/4
                    else:
                        print("Found unexpected base!")
                elif self[x] == 'T':
                    if seq[x] == 'T':
                        distance += 0
                    elif seq[x] in ['A', 'C', 'G', 'R', 'S', 'M', 'V']:
                        distance += 1
                    elif seq[x] in ['Y', 'W', 'K']:
                        distance += 1/2
                    elif seq[x] in ['B', 'D', 'H']:
                        distance += 2/3
                    elif seq[x] in ['N', '-', '?']:
                        distance += 3/4
                    else:
                        print("Found unexpected base!")
                elif self[x] == 'R':
                    if seq[x] in ['C', 'T', 'Y']:
                        distance += 1
                    elif seq[x] in ['A', 'G', 'R']:
                        distance += 1/2
                    elif seq[x] in ['S', 'W', 'K', 'M']:
                        distance += 3/4
                    elif seq[x] in ['D', 'V']:
                        distance += 1/3
                    elif seq[x] in ['B', 'H']:
                        distance += 5/6
                    elif seq[x] in ['N', '-', '?']:
                        distance += 3/4
                    else:
                        print("Found unexpected base!")
                elif self[x] == 'Y':
                    if seq[x] in ['A', 'G', 'R']:
                        distance += 1
                    elif seq[x] in ['C', 'T', 'Y']:
                        distance += 1/2
                    elif seq[x] in ['S', 'W', 'K', 'M']:
                        distance += 3/4
                    elif seq[x] in ['B', 'H']:
                        distance += 2/3
                    elif seq[x] in ['D', 'V']:
                        distance += 5/6
                    elif seq[x] in ['N', '-', '?']:
                        distance += 3/4
                    else:
                        print("Found unexpected base!")
                elif self[x] == 'S':
                    if seq[x] in ['A', 'T', 'W']:
                        distance += 1
                    elif seq[x] in ['C', 'G', 'S']:
                        distance += 1/2
                    elif seq[x] in ['R', 'Y', 'K', 'M']:
                        distance += 3/4
                    elif seq[x] in ['B', 'V']:
                        distance += 2/3
                    elif seq[x] in ['D', 'H']:
                        distance += 5/6
                    elif seq[x] in ['N', '-', '?']:
                        distance += 3/4
                    else:
                        print("Found unexpected base!")
                elif self[x] == 'W':
                    if seq[x] in ['C', 'G', 'S']:
                        distance += 1
                    elif seq[x] in ['A', 'T', 'W']:
                        distance += 1/2
                    elif seq[x] in ['R', 'Y', 'K', 'M']:
                        distance += 3/4
                    elif seq[x] in ['D', 'H']:
                        distance += 2/3
                    elif seq[x] in ['B', 'V']:
                        distance += 5/6
                    elif seq[x] in ['N', '-', '?']:
                        distance += 3/4
                    else:
                        print("Found unexpected base!")
                elif self[x] == 'K':
                    if seq[x] in ['A', 'C', 'M']:
                        distance += 1
                    elif seq[x] in ['G', 'T', 'K']:
                        distance += 1/2
                    elif seq[x] in ['R', 'Y', 'S', 'W']:
                        distance += 3/4
                    elif seq[x] in ['B', 'D']:
                        distance += 2/3
                    elif seq[x] in ['H', 'V']:
                        distance += 5/6
                    elif seq[x] in ['N', '-', '?']:
                        distance += 3/4
                    else:
                        print("Found unexpected base!")
                elif self[x] == 'M':
                    if seq[x] in ['G', 'T', 'K']:
                        distance += 1
                    elif seq[x] in ['A', 'C', 'M']:
                        distance += 1/2
                    elif seq[x] in ['R', 'Y', 'S', 'W']:
                        distance += 3/4
                    elif seq[x] in ['H', 'V']:
                        distance += 2/3
                    elif seq[x] in ['B', 'D']:
                        distance += 5/6
                    elif seq[x] in ['N', '-', '?']:
                        distance += 3/4
                    else:
                        print("Found unexpected base!")
                elif self[x] == 'B':
                    if seq[x] == 'A':
                        distance += 1
                    elif seq[x] in ['C', 'G', 'T', 'Y', 'S', 'K', 'B']:
                        distance += 2/3
                    elif seq[x] in ['R', 'W', 'M']:
                        distance += 5/6
                    elif seq[x] in ['D', 'H', 'V']:
                        distance += 7/9
                    elif seq[x] in ['N', '-', '?']:
                        distance += 3/4
                    else:
                        print("Found unexpected base!")
                elif self[x] == 'D':
                    if seq[x] == 'C':
                        distance += 1
                    elif seq[x] in ['A', 'G', 'T', 'R', 'W', 'K', 'D']:
                        distance += 2/3
                    elif seq[x] in ['Y', 'S', 'M']:
                        distance += 5/6
                    elif seq[x] in ['B', 'H', 'V']:
                        distance += 7/9
                    elif seq[x] in ['N', '-', '?']:
                        distance += 3/4
                    else:
                        print("Found unexpected base!")
                elif self[x] == 'H':
                    if seq[x] == 'G':
                        distance += 1
                    elif seq[x] in ['A', 'C', 'T', 'Y', 'W', 'M', 'H']:
                        distance += 2/3
                    elif seq[x] in ['R', 'S', 'K']:
                        distance += 5/6
                    elif seq[x] in ['B', 'D', 'V']:
                        distance += 7/9
                    elif seq[x] in ['N', '-', '?']:
                        distance += 3/4
                    else:
                        print("Found unexpected base!")
                elif self[x] == 'V':
                    if seq[x] == 'T':
                        distance += 1
                    elif seq[x] in ['A', 'C', 'G', 'R', 'S', 'M', 'H']:
                        distance += 2/3
                    elif seq[x] in ['Y', 'W', 'K']:
                        distance += 5/6
                    elif seq[x] in ['B', 'D', 'H']:
                        distance += 7/9
                    elif seq[x] in ['N', '-', '?']:
                        distance += 3/4
                    else:
                        print("Found unexpected base!")
                else:
                    print("Found unexpected base!")
        return distance


# Expand the class MultipleSeqAlignment in order to add more statistics.
class XMultipleSeqAlignment(MultipleSeqAlignment):

    def get_number_of_records(self, pop=None):
        number_of_records_for_this_pop = 0
        for y in range(0,len(self)):
            if pop == None or pop in self[y].id:
                number_of_records_for_this_pop += 1
        return number_of_records_for_this_pop

    def set_is_haploid(self, haploid):
        self.is_haploid = haploid

    def get_is_haploid(self):
        return self.is_haploid

    def all_records_map_uniquely_to_pops(self, pops=None):
        if pops == None:
            return False
        else:
            unique_map = True
            for record in self:
                count = 0
                for pop in pops:
                    if pop in record.id:
                        count += 1
                if count != 1:
                    unique_map = False
            return unique_map

    def get_number_of_variable_sites(self, pops=[]):
        # Get all sequences of this population.
        seqs = []
        for y in range(0,len(self)):
            if pops == []:
                seqs.append(str(self[y].seq))
            else:
                for pop in pops:
                    if pop in self[y].id:
                        seqs.append(str(self[y].seq))
                        break
        number_of_variable_sites = 0
        for x in range(0,self.get_alignment_length()):
            bases = []
            for y in range(0,len(seqs)):
                if seqs[y][x] in ['a','A','c','C','g','G','t','T']:
                    bases.append(seqs[y][x])
            if len(set(bases)) > 1:
                number_of_variable_sites += 1
        return number_of_variable_sites

    def get_number_of_unique_genotypes(self, pop=None):
        seqs = []
        for record in self:
            if pop == None or pop in record.id:
                seqs.append(str(record.seq))
        return len(set(seqs))

    def get_proportion_of_variable_sites(self, pops=[]):
        return self.get_number_of_variable_sites(pops)/self.get_alignment_length()

    def get_number_of_invariable_sites(self, pops=[]):
        return self.get_alignment_length()-self.get_number_of_variable_sites(pops)

    def get_allele_frequencies(self, x, pops=[]):
        if isinstance(pops, list) == False:
            print("ERROR: Populations are not given as a list (get_allele_frequencies):")
            print(pops)
            sys.exit(1)
        else:
            allele_frequencies = [0, 0, 0, 0]
            for y in range(0,len(self)):
                if pops == []:
                    if self[y].seq[x] is 'A':
                        allele_frequencies[0] += 1
                    elif self[y].seq[x] is 'C':
                        allele_frequencies[1] += 1
                    elif self[y].seq[x] is 'G':
                        allele_frequencies[2] += 1
                    elif self[y].seq[x] is 'T':
                        allele_frequencies[3] += 1
                else:
                    for pop in pops:
                        if pop in self[y].id:
                            if self[y].seq[x] is 'A':
                                allele_frequencies[0] += 1
                            elif self[y].seq[x] is 'C':
                                allele_frequencies[1] += 1
                            elif self[y].seq[x] is 'G':
                                allele_frequencies[2] += 1
                            elif self[y].seq[x] is 'T':
                                allele_frequencies[3] += 1
            return allele_frequencies

    def get_is_biallelic_per_site(self, x, pops=[]):
        if isinstance(pops, list) == False:
            print("ERROR: Populations are not given as a list (get_is_biallelic_per_site):")
            print(pops)
            sys.exit(1)
        else:
            allele_frequencies = [0, 0, 0, 0]
            for pop in pops:
                pop_allele_frequencies = self.get_allele_frequencies(x, [pop])
                allele_frequencies[0] += pop_allele_frequencies[0]
                allele_frequencies[1] += pop_allele_frequencies[1]
                allele_frequencies[2] += pop_allele_frequencies[2]
                allele_frequencies[3] += pop_allele_frequencies[3]
            allele_frequency_zeros = 0
            if allele_frequencies[0] == 0:
                allele_frequency_zeros += 1
            if allele_frequencies[1] == 0:
                allele_frequency_zeros += 1
            if allele_frequencies[2] == 0:
                allele_frequency_zeros += 1
            if allele_frequencies[3] == 0:
                allele_frequency_zeros += 1
            if allele_frequency_zeros == 2:
                return True
            else:
                return False

    def get_is_variable_per_site(self, x, pops=[]):
        if isinstance(pops, list) == False:
            print("ERROR: Populations are not given as a list (get_is_variable_per_site):")
            print(pops)
            sys.exit(1)
        else:
            allele_frequencies = [0, 0, 0, 0]
            for pop in pops:
                pop_allele_frequencies = self.get_allele_frequencies(x, [pop])
                allele_frequencies[0] += pop_allele_frequencies[0]
                allele_frequencies[1] += pop_allele_frequencies[1]
                allele_frequencies[2] += pop_allele_frequencies[2]
                allele_frequencies[3] += pop_allele_frequencies[3]
            allele_frequency_zeros = 0
            if allele_frequencies[0] == 0:
                allele_frequency_zeros += 1
            if allele_frequencies[1] == 0:
                allele_frequency_zeros += 1
            if allele_frequencies[2] == 0:
                allele_frequency_zeros += 1
            if allele_frequencies[3] == 0:
                allele_frequency_zeros += 1
            if allele_frequency_zeros == 3:
                return False
            else:
                return True

    def get_alleles_per_site(self, x, pops=[]):
        if isinstance(pops, list) == False:
            print("ERROR: Populations are not given as a list (get_alleles_per_site):")
            print(pops)
            sys.exit(1)
        else:
            alleles = []
            for y in range(0,len(self)):
                if pops == []:
                    alleles << self[y].seq[x]
                else:
                    for pop in pops:
                        if pop in self[y].id:
                            alleles.append(self[y].seq[x])
            return alleles

    def get_pi_per_site(self, x, pops=[]):
        # pi is the probability that two randomly chosen sequences from the
        # sample have different alleles at a site x.
        # Following Ruegg et al. (2014, Mol Ecol, A role for migration-linked genes and genomic
        # islands in divergence of a songbird), only biallelic (or in this case monomorphic)
        # SNPs are allowed.
        if isinstance(pops, list) == False:
            print("ERROR: Populations are not given as a list (get_pi_per_site):")
            print(pops)
            sys.exit(1)
        else:
            all_allele_frequencies = self.get_allele_frequencies(x, pops)
            two_allele_frequencies = []
            for all_allele_frequency in all_allele_frequencies:
                if all_allele_frequency is not 0:
                    two_allele_frequencies.append(all_allele_frequency)
            if len(two_allele_frequencies) == 0:
                return 0
            elif len(two_allele_frequencies) == 1:
                return 0
            elif len(two_allele_frequencies) == 2:
                # Use r and a as in Ruegg et al. (2014).
                r = two_allele_frequencies[0]
                a = two_allele_frequencies[1]
                numerator = r * a
                denominator = scipy.special.binom((r+a),2)
                pi = numerator/denominator
                return pi
            elif len(two_allele_frequencies) > 2:
                return None

    def get_pi(self, pops=[]):
        pi_per_site_values = []
        l_k = self.get_alignment_length()
        for x in range(self.get_alignment_length()):
            pi_per_site = self.get_pi_per_site(x, pops)
            if pi_per_site is not None:
                pi_per_site_values.append(pi_per_site)
        if l_k == 0:
            return None
        else:
            return sum(pi_per_site_values)/l_k

    def get_F_st(self, pops=[]):
        if len(pops) != 2:
            print("ERROR: Exactly two populations must be specified to calculate pairwise F_st!")
            sys.exit(1)
        else:
            # Initiate six lists that will be needed for Fst calculation.
            # ninds_dup_1 and ninds_dup_2 and the numbers of individuals for pop1 and pop2, per locus.
            ninds_dup_1 = []
            ninds_dup_2 = []
            # p1 and p2 are the frequencies of allele A for pop1 and pop2, per locus.
            p1 = []
            p2 = []
            # oh1 and oh2 are the frequencies of heterozygotes for pop1 and pop2, per locus.
            oh1 = []
            oh2 = []

            # The number of populations, r is set to two.
            r = 2

            # Consider only biallelic loci, and only individuals without missing data.
            for x in range(self.get_alignment_length()):
                if self.get_is_biallelic_per_site(x, pops):

                    # Get all alleles at this site.
                    alleles_pop1 = self.get_alleles_per_site(x, [pops[0]])
                    alleles_pop2 = self.get_alleles_per_site(x, [pops[1]])
                    alleles_both_pops = alleles_pop1 + alleles_pop2
                    unique_alleles_both_pops = sorted(list(set(alleles_both_pops)))
                    good_unique_alleles_both_pops = []
                    for unique_allele_both_pops in unique_alleles_both_pops:
                        if unique_allele_both_pops in ["A", "C", "G", "T"]:
                            good_unique_alleles_both_pops.append(unique_allele_both_pops)

                    # If sequences are haploid, 'diploidize' them by counting each twice.
                    if self.get_is_haploid() == True:
                        tmp = []
                        for allele in alleles_pop1:
                            tmp.append(allele)
                            tmp.append(allele)
                        alleles_pop1 = tmp
                        tmp = []
                        for allele in alleles_pop2:
                            tmp.append(allele)
                            tmp.append(allele)
                        alleles_pop2 = tmp

                    # Get the good (= no missing data) pairs of alleles for both pops, for this site.
                    good_pairs_of_alleles_pop1 = []
                    good_pairs_of_alleles_pop2 = []
                    for z in range(int(len(alleles_pop1)/2)):
                        if alleles_pop1[z*2] in good_unique_alleles_both_pops and alleles_pop1[(z*2)+1] in good_unique_alleles_both_pops:
                            good_pairs_of_alleles_pop1.append([alleles_pop1[z*2],alleles_pop1[(z*2)+1]])
                    for z in range(int(len(alleles_pop2)/2)):
                        if alleles_pop2[z*2] in good_unique_alleles_both_pops and alleles_pop2[(z*2)+1] in good_unique_alleles_both_pops:
                            good_pairs_of_alleles_pop2.append([alleles_pop2[z*2],alleles_pop2[(z*2)+1]])

                    # Get the number of good individuals in both pops based on the number of good pairs of alleles, for this site.
                    ninds_dup_1_this_site = len(good_pairs_of_alleles_pop1)
                    ninds_dup_2_this_site = len(good_pairs_of_alleles_pop2)

                    # Continue only if good individuals are found in both populations, for this site.
                    if ninds_dup_1_this_site > 0 and ninds_dup_2_this_site > 0:

                        # Get the frequencies of the two alleles, for this site.
                        allele_A_occurrences_pop1 = 0
                        allele_A_occurrences_pop2 = 0
                        for good_pair_of_alleles_pop1 in good_pairs_of_alleles_pop1:
                            for good_allele_pop1 in good_pair_of_alleles_pop1:
                                if good_allele_pop1 == good_unique_alleles_both_pops[0]:
                                    allele_A_occurrences_pop1 += 1
                                elif good_allele_pop1 != good_unique_alleles_both_pops[1]:
                                    print("ERROR: Unexpected allele found at site " + str(x) + ": " + good_allele_pop1 + "!")
                                    sys.exit(1)
                        for good_pair_of_alleles_pop2 in good_pairs_of_alleles_pop2:
                            for good_allele_pop2 in good_pair_of_alleles_pop2:
                                if good_allele_pop2 == good_unique_alleles_both_pops[0]:
                                    allele_A_occurrences_pop2 += 1
                                elif good_allele_pop2 != good_unique_alleles_both_pops[1]:
                                    print("ERROR: Unexpected allele found at site " + str(x) + ": " + good_allele_pop2 + "!")
                                    sys.exit(1)
                        p1_this_site = allele_A_occurrences_pop1/(len(good_pairs_of_alleles_pop1)*2)
                        p2_this_site = allele_A_occurrences_pop2/(len(good_pairs_of_alleles_pop2)*2)

                        # Get the frequencies of heterozygotes, for this site.
                        heterozygote_occurrences_pop1 = 0
                        heterozygote_occurrences_pop2 = 0
                        for good_pair_of_alleles_pop1 in good_pairs_of_alleles_pop1:
                            if good_pair_of_alleles_pop1[0] != good_pair_of_alleles_pop1[1]:
                                heterozygote_occurrences_pop1 += 1
                        for good_pair_of_alleles_pop2 in good_pairs_of_alleles_pop2:
                            if good_pair_of_alleles_pop2[0] != good_pair_of_alleles_pop2[1]:
                                heterozygote_occurrences_pop2 += 1
                        oh1_this_site = heterozygote_occurrences_pop1/len(good_pairs_of_alleles_pop1)
                        oh2_this_site = heterozygote_occurrences_pop2/len(good_pairs_of_alleles_pop2)

                        ninds_dup_1.append(ninds_dup_1_this_site)
                        ninds_dup_2.append(ninds_dup_2_this_site)
                        p1.append(p1_this_site)
                        p2.append(p2_this_site)
                        oh1.append(oh1_this_site)
                        oh2.append(oh2_this_site)

            # Calculation of Fst according to Weir & Cockerham (1984), as in method stamppFst of R package StAMPP.
            # n_bar is the average number of individuals in a population, for each locus.
            n_bar = []
            for x in range(len(ninds_dup_1)):
                n_bar.append((ninds_dup_1[x] + ninds_dup_2[x])/r)
            nc = []
            for x in range(len(ninds_dup_1)):
                nc.append((r * n_bar[x]) - (((ninds_dup_1[x]**2) + (ninds_dup_2[x]**2))/(r * n_bar[x])))
            # p_bar is the average sample frequency of allele A in a population, for each locus.
            p_bar = []
            for x in range(len(ninds_dup_1)):
                p_bar.append(((ninds_dup_1[x] * p1[x])/(r * n_bar[x])) + ((ninds_dup_2[x] * p2[x])/(r * n_bar[x])))
            # s_square is the sample variance of allele A frequencies over populations.
            s_square = []
            for x in range(len(ninds_dup_1)):
                s_square.append(((ninds_dup_1[x] * ((p1[x] - p_bar[x])**2))/n_bar[x]) + ((ninds_dup_2[x] * ((p2[x] - p_bar[x])**2))/n_bar[x]))
            # h_bar is the average heterozygote frequency for allele A.
            h_bar = []
            for x in range(len(ninds_dup_1)):
                h_bar.append(((ninds_dup_1[x] * oh1[x])/(r * n_bar[x])) + ((ninds_dup_2[x] * oh2[x])/(r * n_bar[x])))
            # Equation 2 in WC84.
            a = []
            for x in range(len(ninds_dup_1)):
                if n_bar[x] > 1:
                    a_cand = (n_bar[x]/nc[x]) * (s_square[x] - (1/(n_bar[x] - 1)) * ((p_bar[x] * (1 - p_bar[x])) - (((r - 1)/r) * s_square[x]) - ((1/4) * h_bar[x])))
                    if not math.isnan(a_cand):
                        a.append(a_cand)
            # Equation 3 in WC84.
            b = []
            for x in range(len(ninds_dup_1)):
                if n_bar[x] > 1:
                    b_cand = (n_bar[x]/(n_bar[x] - 1)) * ((p_bar[x] * (1 - p_bar[x])) - (((r - 1)/r) * s_square[x]) - (((2 * n_bar[x] - 1)/(4 * n_bar[x])) * h_bar[x]))
                    if not math.isnan(b_cand):
                        b.append(b_cand)
            # Equation 4 in WC84.
            c = []
            for x in range(len(ninds_dup_1)):
                if n_bar[x] > 1:
                    c_cand = (1/2) * h_bar[x]
                    if not math.isnan(c_cand):
                        c.append(c_cand)
            if (sum(a) + sum(b) + sum(c)) == 0:
                return None
            else:
                fst = sum(a)/(sum(a) + sum(b) + sum(c))
                return fst

    def get_d_xy(self, pops=[]):
        if len(pops) != 2:
            print("ERROR: Exactly two populations must be specified to calculate pairwise d_xy!")
            sys.exit(1)
        else:
            l_k = self.get_alignment_length()
            sum_of_quotients = 0
            for x in range(self.get_alignment_length()):
                if self.get_is_biallelic_per_site(x, pops):
                    all_allele_frequencies0 = self.get_allele_frequencies(x, [pops[0]])
                    all_allele_frequencies1 = self.get_allele_frequencies(x, [pops[1]])
                    if sum(all_allele_frequencies0) > 0 and sum(all_allele_frequencies1) > 0:
                        two_allele_frequencies0 = []
                        two_allele_frequencies1 = []
                        if all_allele_frequencies0[0] != 0 or all_allele_frequencies1[0] != 0:
                            two_allele_frequencies0.append(all_allele_frequencies0[0])
                            two_allele_frequencies1.append(all_allele_frequencies1[0])
                        if all_allele_frequencies0[1] != 0 or all_allele_frequencies1[1] != 0:
                            two_allele_frequencies0.append(all_allele_frequencies0[1])
                            two_allele_frequencies1.append(all_allele_frequencies1[1])
                        if all_allele_frequencies0[2] != 0 or all_allele_frequencies1[2] != 0:
                            two_allele_frequencies0.append(all_allele_frequencies0[2])
                            two_allele_frequencies1.append(all_allele_frequencies1[2])
                        if all_allele_frequencies0[3] != 0 or all_allele_frequencies1[3] != 0:
                            two_allele_frequencies0.append(all_allele_frequencies0[3])
                            two_allele_frequencies1.append(all_allele_frequencies1[3])
                        if len(two_allele_frequencies0) != 2:
                            print("ERROR: Wrong number of allele frequencies!")
                            print(all_allele_frequencies0)
                            print(all_allele_frequencies1)
                            sys.exit(1)
                        r0 = two_allele_frequencies0[0]
                        a0 = two_allele_frequencies0[1]
                        r1 = two_allele_frequencies1[0]
                        a1 = two_allele_frequencies1[1]
                        numerator = r0*a1 + r1*a0
                        denominator = (r0+a0) * (r1+a1)
                        sum_of_quotients += numerator/denominator
            d_xy = sum_of_quotients/l_k
            return d_xy

    def get_d_f(self, pops=[]):
        if len(pops) != 2:
            print("ERROR: Exactly two populations must be specified to calculate pairwise d_f!")
            sys.exit(1)
        else:
            l_k = self.get_alignment_length()
            number_of_fixed_snps = 0
            for x in range(self.get_alignment_length()):
                if self.get_is_biallelic_per_site(x, pops):
                    pop0_variable = self.get_is_variable_per_site(x, [pops[0]])
                    pop1_variable = self.get_is_variable_per_site(x, [pops[1]])
                    if pop0_variable == False and pop1_variable == False:
                        number_of_fixed_snps += 1
            d_f = number_of_fixed_snps/l_k
            return d_f


# Parse the command line arguments.
parser = argparse.ArgumentParser(
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description=textwrap.dedent('''\
      %(prog)s
    -----------------------------------------
      Reads nexus formatted strings and produces
      haplotype genealogy graphs using the Fitch algorithm.
      Start e.g. with
      fitchi.py example.nex example.html -p pop3 pop5
      Info: www.evoinformatics.eu/fitchi.htm
    '''))
parser.add_argument(
    '-v', '--version',
    action='version',
    version='%(prog)s 1.1')
parser.add_argument(
    '-p', '--populations',
    nargs='*',
    type=str,
    help="One or more population identifiers (default: none).")
parser.add_argument(
    '-f', '--from',
    nargs=1,
    type=int,
    default=[1],
    dest='start',
    help="Start position of analysis window (count starts at 1) (default: 1)."
    )
parser.add_argument(
    '-t', '--to',
    nargs=1,
    type=int,
    default=[-1],
    dest='end',
    help="End position of analysis window (count starts at 1) (default: last position of alignment)."
    )
parser.add_argument(
    '-e', '--min-edge-length',
    nargs=1,
    type=int,
    default=[1],
    dest='min_edge_length',
    help="Minimum edge length for display in haplotype genealogy graph (default: 1)."
    )
parser.add_argument(
    '-n', '--min-node-size',
    nargs=1,
    type=int,
    default=[1],
    dest='min_node_size',
    help="Minimum node size for display in haplotype genealogy graph (default: 1)."
    )
parser.add_argument(
    '-x', '--transversions-only',
    action='store_true',
    dest='transversions_only',
    help="Ignore transitions and show transversions only (default: off)."
    )
parser.add_argument(
    '--haploid',
    action='store_true',
    dest='haploid',
    help="Sequences are haploid (default: off). This only affects Fst calculations. If not specified, each pair of two consecutive sequences is assumed to be from the same individual."
    )
parser.add_argument(
    '-m', '--radius-multiplier',
    nargs=1,
    type=float,
    default=[1.0],
    dest='radius_multiplier',
    help="Scale factor for the size of node radi (default: 1.0).")
parser.add_argument(
    '-s', '--seed',
    nargs=1,
    type=int,
    default=[-1],
    dest='seed',
    help="Specifies a random number seed."
    )
parser.add_argument(
    'infile',
    nargs='?',
    type=argparse.FileType('r'),
    default='-',
    help='The input file name.')
parser.add_argument(
    'outfile', nargs='?',
    type=argparse.FileType('w'),
    default=sys.stdout,
    help='The output file name.')
args = parser.parse_args()
infile = args.infile
outfile = args.outfile
window_start_pos = args.start[0]-1
window_end_pos = args.end[-1]
minimum_edge_length = args.min_edge_length[0]
minimum_node_size = args.min_node_size[0]
radius_multiplier = args.radius_multiplier[0]
seed = args.seed[0]
transversions_only = args.transversions_only
haploid = args.haploid

# Initialize the random number generator if a seed value has been provided.
if seed == -1:
    seed = random.randint(0, 99999)
random.seed(seed)

# Make sure sensible values are specified for the window start and end.
if window_start_pos < 0:
    print("ERROR: The start position of the analysis window must be at least 1!")
    sys.exit(1)
elif window_end_pos != -1:
    if window_end_pos <= window_start_pos:
        print("ERROR: The end position of the analysis window must be greater than the start position!")
        sys.exit(1)
pops = args.populations
if infile.isatty():
    print("No input file specified, and no input piped through stdin!")
    print("Use '-h' to see available options.")
    sys.exit(1)

# Define a color scheme.
# Colors use the Solarized color scheme of http://ethanschoonover.com/solarized.
if pops == None:
    pops = []
colors = []
if len(pops) == 1:
    # base01
    colors = ['586e75']
elif len(pops) == 2:
    # red, cyan
    colors = ['dc322f', '2aa198']
elif len(pops) == 3:
    # red, cyan, violet
    colors = ['dc322f', '2aa198', '6c71c4']
elif len(pops) == 4:
    # red, cyan, violet, yellow
    colors = ['dc322f', '2aa198', '6c71c4', 'b58900']
elif len(pops) == 5:
    # red, cyan, violet, yellow, green
    colors = ['dc322f', '2aa198', '6c71c4', 'b58900', '859900']
elif len(pops) == 6:
    # red, cyan, violet, yellow, green, magenta
    colors = ['dc322f', '2aa198', '6c71c4', 'b58900', '859900', 'd33682']
elif len(pops) == 7:
    # red, cyan, violet, yellow, green, magenta, blue
    colors = ['dc322f', '2aa198', '6c71c4', 'b58900', '859900', 'd33682', '268bd2']
elif len(pops) == 8:
    # yellow, green, cyan, blue, violet, magenta, red, orange
    colors = ['859900', 'b58900', '2aa198', '268bd2', '6c71c4', 'd33682', 'dc322f', 'cb4b16']
elif len(pops) == 9:
    # yellow, green, cyan, blue, violet, magenta, red, orange, base03
    colors = ['859900', 'b58900', '2aa198', '268bd2', '6c71c4', 'd33682', 'dc322f', 'cb4b16', '002b36']
elif len(pops) > 9:
    # yellow, green, cyan, blue, violet, magenta, red, orange, base03, base01
    colors = ['859900', 'b58900', '2aa198', '268bd2', '6c71c4', 'd33682', 'dc322f', 'cb4b16', '002b36', '586e75']
# base1
rest_color = '93a1a1'

# Parse the input.
align = None
tree = None
inlines = infile.readlines()
if inlines[0][0:6].lower() == '#nexus':
    # Assume the input is in nexus format. Maximally one tree string is read.
    in_matrix = False
    in_tree = False
    records = []
    for line in inlines:
        clean_line = line.strip()
        if "[" in clean_line and "]" in clean_line:
            tmp = ""
            in_comment = False
            for letter in clean_line:
                if letter == "[":
                    in_comment = True
                elif letter == "]":
                    in_comment = False
                elif in_comment == False:
                    tmp += letter
            clean_line = tmp
        if clean_line.lower() == 'matrix':
            in_matrix = True
        elif clean_line == ';':
            in_matrix = False
            in_tree = False
        elif "format" in clean_line.lower():
            if "interleave" in clean_line.lower():
                print("ERROR: Could not parse the alignment (should be sequential nexus format, but the format specification says that it is interleaved)!")
                sys.exit(1)
        elif in_matrix and clean_line is not '':
            line_ary = clean_line.split()
            if len(line_ary) == 1:
                print("ERROR: Could not parse the alignment (should be sequential nexus format, but looks like it is interleaved)!")
                sys.exit(1)
            elif len(line_ary) > 2:
                print("ERROR: Could not parse the alignment!")
                sys.exit(1)
            else:
                seq_string = line_ary[1].upper()
            pattern = re.compile("^[a-zA-Z0-9_\.\-]+?$")
            hit = pattern.search(line_ary[0])
            if hit == None:
                print("ERROR: Taxon labels should include only 'A'-'Z', 'a'-'z', '0'-'9', '.', _', and '-'! Offending taxon label: " + line_ary[0] + ".")
                sys.exit(1)
            if window_end_pos == -1:
                seq_string = seq_string[window_start_pos:]
            else:
                seq_string = seq_string[window_start_pos:window_end_pos]
            records.append(
                SeqRecord(
                    Seq(seq_string,
                        generic_dna),
                        id = line_ary[0]))
        elif line.strip() == 'begin trees;':
            in_tree = True
        elif line.strip() == 'end;':
            in_tree = False
        elif in_tree and line.strip() is not '':
            tree_string_raw = line
            tree_patterns = re.search('\(.+\)',tree_string_raw)
            tree_string = tree_patterns.group(0)
            tree = Tree(tree_string)
    if records == []:
        print("ERROR: File could not be parsed!")
    else:
        for record in records:
            if record.id not in tree_string:
                print("ERROR: Record id " + record.id + " not found in tree string!")
                sys.exit(0)

    align = XMultipleSeqAlignment(records)
    if haploid:
        align.set_is_haploid(True)
    else:
        align.set_is_haploid(False)
else:
    print("ERROR: Unexpected file format!")
    sys.exit(1)

# Parse the newick tree string.
tree.parse_newick_string(pops)

# Assign sequences to terminal nodes.
nodes = tree.get_nodes()
for node in nodes:
    for seq in align:
        if node.get_id() == seq.id:
            node.set_sequences([str(seq.seq.upper())])
            break

# Make sure all non-internal nodes have sequences.
all_seqs_found = True
for node in nodes:
    node_id = node.get_id()
    if node_id[:12] != 'internalNode':
        if node.get_sequences() == []:
            print("ERROR: No sequence was found for node " + node_id + "!")
            all_seqs_found = False
if all_seqs_found == False:
    sys.exit(1)

# Reconstruct ancestral sequences using the Fitch algorithm.
tree.reconstruct_ancestral_sequences()

# Calculate Fitch distances.
tree.calculate_fitch_distances(transversions_only)

# Find the extant progeny for each edge.
tree.assign_progeny_ids()

# Calculate gsi values before the tree is reduced.
gsis = []
for pop in pops:
    gsis.append(tree.get_gsi(pop))

# Reduce the tree.
tree.reduce(minimum_edge_length, minimum_node_size)

# Position the tree.
tree.position('neato', minimum_node_size, radius_multiplier)

# Produce the svg tree string.
svg_string = tree.to_svg(840, 700, 10, minimum_node_size, radius_multiplier, colors, rest_color)

# Initiate the html output string.
html_string = ''
html_string += '<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"\n'
html_string += '"http://www.w3.org/TR/html4/loose.dtd">\n'
html_string += '<html>\n'
html_string += '  <head>\n'
html_string += '    <title>Fitchi results</title>\n'
html_string += '    <meta http-equiv="Content-Type" content="text/html; charset=iso-8859-1">\n'
html_string += '\n'
html_string += '    <style type="text/css">\n'
html_string += '      a:link { text-decoration:none; color:#000000; }\n'
html_string += '      a:visited { text-decoration:none; color:#000000; }\n'
html_string += '      a:hover { text-decoration:none; color:#ffffff; background-color:#000000; }\n'
html_string += '      a:active { text-decoration:none; color:#ffffff; background-color:#000000; }\n'
html_string += '      a:focus { text-decoration:none; color:#ffffff; background-color:#000000; }\n'
html_string += '      #relativeSVG { position:relative; width:840px; height:236px; z-index:2 }\n'
html_string += '      #absoluteAst { position:absolute; top:20px; left:20px; width:20px; height:20px; z-index:1; background-color:#e05030 }\n'
html_string += '      td { font-family:helvetica; font-size:12px }\n'
html_string += '      tr.spaceUnder > td { padding-bottom: 1em; }\n'
html_string += '      tr.doubleSpaceUnder > td { padding-bottom: 2em; }\n'
html_string += '      tr.largeSpaceUnder > td { padding-bottom: 8em; }\n'
html_string += '      tr.smallSpaceUnder > td { padding-bottom: 0.2em; }\n'
html_string += '      tr.spaceOver > td { padding-top: 1em; }\n'
html_string += '      tr.spaceOverAndLargeSpaceUnder > td { padding-top: 1em; padding-bottom: 8em; }\n'
html_string += '    </style>\n'
html_string += '    <script type="text/javascript">\n'
html_string += '      <!--\n'
html_string += '        function legend_toggle () {\n'
html_string += '          if(document.getElementById("legend").style.display == "none") {\n'
html_string += '            document.getElementById("legend").style.display = "inline";\n'
html_string += '          } else {\n'
html_string += '            document.getElementById("legend").style.display = "none";\n'
html_string += '          }\n'
html_string += '          if(document.getElementById("show").style.display == "none") {\n'
html_string += '            document.getElementById("show").style.display = "inline";\n'
html_string += '          } else {\n'
html_string += '            document.getElementById("show").style.display = "none";\n'
html_string += '          }\n'
html_string += '          if(document.getElementById("hide").style.display == "none") {\n'
html_string += '            document.getElementById("hide").style.display = "inline";\n'
html_string += '          } else {\n'
html_string += '            document.getElementById("hide").style.display = "none";\n'
html_string += '          }\n'
html_string += '        }\n'
html_string += '      //-->\n'
html_string += '    </script>\n'
html_string += '  </head>\n'
html_string += '\n'
html_string += '  <body>\n'
html_string += '    <div align="center">\n'
html_string += '      <table width="840" border="0" cellpadding="0" cellspacing="0">\n'
html_string += '        <tr>\n'
html_string += '          <td style="font-family:helvetica; font-size:54px; font-weight:bold">\n'
html_string += '            <svg width="50" height="50">\n'
html_string += '              <defs>\n'
html_string += '                <radialGradient id="radgrad" cx=".9" cy=".1" r="1">\n'
html_string += '                  <stop  offset="0" style="stop-color:grey"/>\n'
html_string += '                  <stop  offset="1" style="stop-color:black"/>\n'
html_string += '                </radialGradient>\n'
html_string += '                <mask id="m_0"><circle fill="url(#radgrad)" cx="30.964" cy="31.513" r="11.0695"/></mask>\n'
html_string += '                <mask id="m_1"><circle fill="url(#radgrad)" cx="39.693" cy="7.257" r="6.257"/></mask>\n'
html_string += '                <mask id="m_2"><circle fill="url(#radgrad)" cx="43.419" cy="45.559" r="3.44"/></mask>\n'
html_string += '                <mask id="m_3"><circle fill="url(#radgrad)" cx="7.91" cy="39.167" r="4.769"/></mask>\n'
html_string += '              </defs>\n'
html_string += '              <path fill="#94A2A1" stroke="#94A2A1" stroke-width="0.5" d="M32.802,20.593c5.745,0.957,9.882,6.31,9.147,12.25c-0.778,6.312-6.659,10.731-13.029,9.548c-6.298-1.18-10.148-7.31-8.733-13.371C21.521,23.296,27.118,19.643,32.802,20.593z"/>\n'
html_string += '              <line fill="none" stroke="#94A2A1" stroke-width="0.5" x1="39.693" y1="7.257" x2="30.965" y2="30.919"/>\n'
html_string += '              <line fill="none" stroke="#94A2A1" stroke-width="0.5" x1="30.965" y1="30.919" x2="7.91" y2="39.04"/>\n'
html_string += '              <line fill="none" stroke="#94A2A1" stroke-width="0.5" x1="30.965" y1="30.919" x2="43.419" y2="45.559"/>\n'
html_string += '              <path fill="#2BA199" d="M35.684,21.532c4.153,1.969,6.83,6.438,6.229,11.311c-0.242,1.953-1.007,3.819-2.22,5.4l-8.729-7.324L35.684,21.532z"/>\n'
html_string += '              <path fill="#859B3B" d="M30.965,30.919l1.834-10.292c1.02,0.17,1.986,0.479,2.885,0.905L30.965,30.919z"/>\n'
html_string += '              <path fill="#2E8BCB" d="M30.965,30.919l8.729,7.324c-1.652,2.136-4.045,3.601-6.723,4.11L30.965,30.919z"/>\n'
html_string += '              <path fill="#DC342E" d="M20.376,28.451c0.462-1.582,1.273-3.029,2.349-4.244l8.24,6.711L20.376,28.451z"/>\n'
html_string += '              <path fill="#D43883" d="M30.965,30.919l-6.627,9.429c-2.32-1.742-3.945-4.395-4.32-7.505c-0.184-1.522-0.045-3.01,0.359-4.392L30.965,30.919z"/>\n'
html_string += '              <path fill="#6F71B5" d="M28.059,42.163c-1.356-0.367-2.618-0.989-3.721-1.815l6.627-9.429l2.008,11.437C31.268,42.668,29.6,42.58,28.059,42.163z"/>\n'
html_string += '              <path fill="#CB4E27" d="M30.965,30.919l-8.24-6.711c1.634-1.845,3.878-3.152,6.443-3.581c0.66-0.114,1.326-0.164,1.993-0.153L30.965,30.919z"/>\n'
html_string += '              <path fill="#B48B2E" d="M31.161,20.474c0.548,0.008,1.096,0.061,1.638,0.153l-1.834,10.292L31.161,20.474z"/>\n'
html_string += '              <circle fill="#94A2A1" stroke="#94A2A1" stroke-width="0.5" cx="39.693" cy="7.257" r="6.257"/>\n'
html_string += '              <circle fill="#CB4E27" cx="39.693" cy="7.257" r="6.257"/>\n'
html_string += '              <path fill="#6F71B5" d="M39.693,7.257l2.311,5.816c-0.715,0.283-1.494,0.44-2.311,0.44c-0.836,0-1.635-0.165-2.362-0.461L39.693,7.257z"/>\n'
html_string += '              <path fill="#2BA199" d="M39.693,7.257l5.438,3.093c-0.698,1.228-1.803,2.196-3.129,2.723L39.693,7.257z"/>\n'
html_string += '              <path fill="#859B3B" d="M39.693,7.257V1c3.455,0,6.257,2.801,6.257,6.257c0,1.124-0.298,2.18-0.817,3.093L39.693,7.257z"/>\n'
html_string += '              <circle fill="#94A2A1" stroke="#94A2A1" stroke-width="0.5" cx="43.419" cy="45.559" r="3.44"/>\n'
html_string += '              <circle fill="#6F71B5" cx="43.419" cy="45.559" r="3.44"/>\n'
html_string += '              <circle fill="#94A2A1" stroke="#94A2A1" stroke-width="0.5" cx="7.91" cy="39.167" r="4.769"/>\n'
html_string += '              <circle fill="#2BA199" cx="7.91" cy="39.167" r="4.769"/>\n'
html_string += '              <path fill="#DC342E" d="M7.91,39.04l-1.548,4.64c-1.873-0.643-3.221-2.42-3.221-4.513c0-2.634,2.135-4.769,4.769-4.769V39.04z"/>\n'
html_string += '              <circle fill="#94A2A1" stroke="#94A2A1" stroke-width="0.5" cx="16.361" cy="36.069" r="0.255"/>\n'
html_string += '              <circle fill="#94A2A1" stroke="#94A2A1" stroke-width="0.5" cx="36.061" cy="17.115" r="0.255"/>\n'
html_string += '              <circle fill="white" cx="30.964" cy="31.513" r="11.0695" mask="url(#m_0)"/>\n'
html_string += '              <circle fill="white" cx="39.693" cy="7.257" r="6.257" mask="url(#m_1)"/>\n'
html_string += '              <circle fill="white" cx="43.419" cy="45.559" r="3.44" mask="url(#m_2)"/>\n'
html_string += '              <circle fill="white" cx="7.91" cy="39.167" r="4.769" mask="url(#m_3)"/>\n'
html_string += '            </svg>\n'
html_string += '            <a name="Fitchi" href="http://www.evoinformatics.eu" style="color:#000000; text-decoration:none; background-color:#ffffff">Fitchi</a><br><br>\n'
html_string += '          </td>\n'
html_string += '        </tr>\n'

# The summary section.
html_string += '        <tr class="smallSpaceUnder">\n'
html_string += '          <td style="font-size:30px; font-weight:bold">Summary</td>\n'
html_string += '        </tr>\n'
html_string += '        <tr class="largeSpaceUnder">\n'
html_string += '          <td>\n'
html_string += '            <table width="840" border="0" cellpadding="0" cellspacing="1">\n'
html_string += '              <tr>\n'
html_string += '                <td width="160" style="font-weight:bold">File Name</td>\n'
if infile.name == '<stdin>':
    html_string += '                <td>STDIN</td>\n'
else:
    html_string += '                <td>' + str(infile.name) + '</td>\n'
html_string += '              </tr>\n'
html_string += '              <tr>\n'
html_string += '                <td width="160" style="font-weight:bold"># Sites</td>\n'
html_string += '                <td>' + str(align.get_alignment_length())
if window_start_pos != 0 or window_end_pos != -1:
    html_string += ' (positions '
    html_string += str(window_start_pos+1) + '-'
    if window_end_pos == -1:
        html_string += str(align.get_alignment_length()+window_start_pos)
    else:
        html_string += str(window_end_pos)
    html_string += ')'
html_string += '</td>\n'
html_string += '              </tr>\n'
html_string += '              <tr>\n'
html_string += '                <td width="160" style="font-weight:bold"># Sequence records</td>\n'
html_string += '                <td>' + str(align.get_number_of_records())
if pops != []:
    pop_string = ' ('
    for pop in pops:
        pop_string += str(align.get_number_of_records(pop)) + ' x ' + pop + ', '
    pop_string = pop_string[:-2]
    html_string += pop_string + ')'
html_string += '</td>\n'
html_string += '              </tr>\n'
html_string += '            </table>\n'
html_string += '          </td>\n'
html_string += '        </tr>\n'

# The haplotype genealogy section.
html_string += '        <tr class="smallSpaceUnder">\n'
html_string += '          <td style="font-size:30px; font-weight:bold"><a href="http://onlinelibrary.wiley.com/doi/10.1111/j.1365-294X.2011.05066.x/abstract">Haplotype genealogy</a></td>\n'
html_string += '        </tr>\n'
html_string += '        <tr class="spaceUnder">\n'
html_string += '          <td>\n'
html_string += '            <table width="840" border="0" cellpadding="0" cellspacing="1">\n'

html_string += '              <tr>\n'
html_string += '                <td width="160" style="font-weight:bold">Node size</td>\n'
html_string += '                <td># Sequence records</td>\n'
html_string += '              </tr>\n'
html_string += '              <tr>\n'
html_string += '                <td width="160" style="font-weight:bold">Edge length</td>\n'
if transversions_only:
    html_string += '                <td># transversions</td>\n'
else:
    html_string += '                <td># substitutions (transitions or transversions)</td>\n'
html_string += '              </tr>\n'


html_string += '              <tr>\n'
html_string += '                <td width="160" style="font-weight:bold">Minimum node size</td>\n'
html_string += '                <td>' + str(minimum_node_size) + ' sequence record'
if minimum_node_size > 1:
    html_string += 's'
html_string += '</td>\n'
html_string += '              </tr>\n'
html_string += '              <tr>\n'
html_string += '                <td width="160" style="font-weight:bold">Minimum edge length</td>\n'
if transversions_only:
    html_string += '                <td>' + str(minimum_edge_length) + ' transversion'
else:
    html_string += '                <td>' + str(minimum_edge_length) + ' substitution'
if minimum_edge_length > 1:
    html_string += 's'
html_string += '</td>\n'
html_string += '              </tr>\n'
html_string += '              <tr>\n'
html_string += '                <td width="160" style="font-weight:bold"># nodes</td>\n'
html_string += '                <td>' + str(tree.get_number_of_nodes()) + '</td>\n'
html_string += '              </tr>\n'
html_string += '              <tr>\n'
html_string += '                <td width="160" style="font-weight:bold"># edges</td>\n'
html_string += '                <td>' + str(tree.get_number_of_edges()) + '</td>\n'
html_string += '              </tr>\n'
html_string += '              <tr>\n'
html_string += '                <td width="160" style="font-weight:bold">Total Fitch distance</td>\n'
html_string += '                <td>'
total_fitch_distance = 0
for edge in tree.get_edges():
    total_fitch_distance += edge.get_fitch_distance()
if transversions_only:
    html_string += str(total_fitch_distance) + ' transversions</td>\n'
else:
    html_string += str(total_fitch_distance) + ' substitutions</td>\n'
html_string += '              </tr>\n'
html_string += '              <tr>\n'
html_string += '                <td width="160" style="font-weight:bold">Random number seed</td>\n'
html_string += '                <td>' + str(seed) + '</td>\n'
html_string += '              </tr>\n'
html_string += '            </table>\n'
html_string += '          </td>\n'
html_string += '        </tr>\n'
html_string += '        <tr>\n'
html_string += '          <td align="center" style="border: 1px solid black">\n'
svg_lines = svg_string.split('\n')
for line in svg_lines:
    html_string += '            ' + line + '\n'
html_string += '          </td>\n'
html_string += '        </tr>\n'
if len(tree.get_nodes()) > 0:
    html_string += '        <tr class="spaceOver">\n'
    html_string += '          <td style="font-weight:bold">\n'
    html_string += '            <span id="show" style="display:inline">\n'
    html_string += '              <button onclick="legend_toggle()">Show legend</button>\n'
    html_string += '            </span>\n'
    html_string += '            <span id="hide" style="display:none">\n'
    html_string += '              <button onclick="legend_toggle()">Hide legend</button>\n'
    html_string += '            </span>\n'
    html_string += '          </td>\n'
    html_string += '        </tr>\n'
html_string += '        <tr class="spaceOverAndLargeSpaceUnder">\n'
html_string += '          <td align="center">\n'
html_string += '            <span id="legend" style="display:none;">\n'
html_string += '              <div style="border-width:1px; border-style:solid; border-color:#000000;">\n'
if pops != []:
    html_string += '                <table width="800" cellpadding="0" cellspacing="1">\n'
    html_string += '                  <tr class="spaceOver">\n'
    html_string += '                    <td width="160" style="font-weight:bold; font-family:Courier">Population</td>\n'
    html_string += '                    <td width="80" style="font-weight:bold; font-family:Courier" colspan="2">Color</td>\n'
    html_string += '                    <td width="250" style="font-weight:bold; font-family:Courier"># Nodes with population presence</td>\n'
    html_string += '                    <td width="210" style="font-weight:bold; font-family:Courier">% Presence in non-empty nodes</td>\n'
    html_string += '                  </tr>\n'
    pop_count = 0
    for pop in pops:
        html_string += '                  <tr>\n'
        html_string += '                    <td width="160" style="font-family:Courier">' + pop + '</td>\n'
        if pop_count >= len(colors):
            html_string += '                    <td width="40" bgcolor="#' + rest_color + '"></td>\n'
        else:
            html_string += '                    <td width="40" bgcolor="#' + colors[pop_count] + '"></td>\n'
        html_string += '                    <td width="40"></td>\n'
        node_count_with_pop = 0
        non_empty_node_count = 0
        for node in tree.get_nodes():
            if pop in node.get_pops():
                node_count_with_pop += 1
            if node.get_size() > 0:
                non_empty_node_count += 1
        html_string += '                    <td style="font-family:Courier">' + str(node_count_with_pop) + '</td>\n'
        if non_empty_node_count > 0:
            html_string += '                    <td style="font-family:Courier">' + "{0:.2f}".format(100*node_count_with_pop/non_empty_node_count) + '</td>\n'
        else:
            html_string += '                    <td style="font-family:Courier">NA</td>\n'
        html_string += '                  </tr>\n'
        pop_count += 1
    html_string += '                  <tr class="doubleSpaceUnder">\n'
    html_string += '                    <td colspan="5"></td>\n'
    html_string += '                  </tr>\n'
    html_string += '                </table>\n'
html_string += '                <table width="800" cellpadding="0" cellspacing="1">\n'
html_string += '                  <!-- Legend: string start -->\n'
nodes = tree.get_nodes()
node_count = 0
first_node = True
for node in nodes:
    if first_node == True:
        first_node = False
        html_string += '                  <tr class="spaceOver">\n'
    else:
        html_string += '                  <tr>\n'
    html_string += '                    <td width="160" style="font-weight:bold; font-family:Courier">Node ' + str(node_count+1) + '</td>\n'
    html_string += '                    <td></td>\n'
    html_string += '                  </tr>\n'
    html_string += '                  <tr>\n'
    html_string += '                    <td valign="top" style="font-family:Courier">Sequence'
    if len(set(node.get_sequences())) > 1:
        html_string += 's'
    html_string += ': </td>\n'
    html_string += '                    <td style="font-family:Courier">\n'
    html_string += '                      <div style="width: 640px; overflow: auto;">\n'
    sequence_string = ''
    for sequence in set(node.get_sequences()):
        sequence_string += '                      ' + sequence + ',<br>\n'
    sequence_string = sequence_string[:-6]
    html_string += sequence_string
    html_string += '\n'
    html_string += '                      </div>\n'
    html_string += '                    </td>\n'
    html_string += '                  </tr>\n'
    html_string += '                  <tr>\n'
    html_string += '                    <td style="font-family:Courier">Size: </td>\n'
    html_string += '                    <td style="font-family:Courier">' + str(node.get_size()) + ' sequence record'
    if node.get_size() > 1:
        html_string += 's'
    if node.get_size() > 0 and pops != []:
        per_pop_sizes = node.get_per_pop_sizes()
        html_string += ' ('
        pops_string = ''
        for x in range(len(pops)):
            if per_pop_sizes[x] > 0:
                pops_string += str(per_pop_sizes[x]) + ' x ' + pops[x] + ', '
        pops_string = pops_string[:-2]
        html_string += pops_string
        html_string += ')'
    html_string += '</td>\n'
    html_string += '                  </tr>\n'
    html_string += '                  <tr class="spaceUnder">\n'
    html_string += '                    <td valign="top" style="font-family:Courier">Sequence record'
    if len(node.get_record_ids()) > 1:
        html_string += 's'
    html_string += ': </td>\n'
    html_string += '                    <td style="font-family:Courier">\n'
    if len(node.get_record_ids()) > 0:
        record_ids_string = ''
        for record_id in sorted(node.get_record_ids()):
            record_ids_string += '                    ' + record_id + ',<br>\n'
        record_ids_string = record_ids_string[:-6]
        html_string += record_ids_string + '\n'
    else:
        html_string += '                      None'
    html_string += '                    </td>\n'
    html_string += '                  </tr>\n'
    node_count += 1
html_string += '                  <!-- Legend: string end -->\n'
html_string += '                </table>\n'
html_string += '              </div>\n'
html_string += '            </span>\n'
html_string += '          </td>\n'
html_string += '        <tr>\n'

# The nucleotide diversity section.
html_string += '        <tr class="smallSpaceUnder">\n'
html_string += '          <td style="font-size:30px; font-weight:bold">Nucleotide diversity</td>\n'
html_string += '        </tr>\n'
html_string += '        <tr class="largeSpaceUnder">\n'
html_string += '          <td style="font-family:helvetica; font-size:12px">\n'
html_string += '            <table width="840" border="0" cellpadding="0" cellspacing="1">\n'
html_string += '              <tr>\n'
html_string += '                <td width="168" style="font-weight:bold">Population</td>\n'
html_string += '                <td width="168" style="font-weight:bold">Variable sites</td>\n'
html_string += '                <td width="168" style="font-weight:bold">Invariable sites</td>\n'
html_string += '                <td width="168" style="font-weight:bold">Proportion variable</td>\n'
html_string += '                <td width="168" style="font-weight:bold">&pi;</td>\n'
html_string += '              </tr>\n'
html_string += '              <tr>\n'
html_string += '                <td width="168" style="font-weight:bold">All</td>\n'
html_string += '                <!-- Total variable: string start -->\n'
html_string += '                <td width="168">' + str(align.get_number_of_variable_sites()) + '</td>\n'
html_string += '                <!-- Total variable: string end -->\n'
html_string += '                <td width="168">' + str(align.get_number_of_invariable_sites()) + '</td>\n'
html_string += '                <!-- Proportion variable: string start -->\n'
html_string += '                <td width="168">' + "{0:.4f}".format(align.get_proportion_of_variable_sites()) + '</td>\n'
html_string += '                <!-- Proportion variable: string end -->\n'
html_string += '                <!-- Pi: string start -->\n'
html_string += '                <td width="168">' + "{0:.4f}".format(align.get_pi()) + '</td>\n'
html_string += '                <!-- Pi: string end -->\n'
html_string += '              </tr>\n'
if pops != None:
    for pop in pops:
        html_string += '              <tr>\n'
        html_string += '                <td width="168" style="font-weight:bold">' + pop + '</td>\n'
        html_string += '                <td width="168">' + str(align.get_number_of_variable_sites([pop])) + '</td>\n'
        html_string += '                <td width="168">' + str(align.get_number_of_invariable_sites([pop])) + '</td>\n'
        html_string += '                <td width="168">' + "{0:.4f}".format(align.get_proportion_of_variable_sites([pop])) + '</td>\n'
        html_string += '                <td width="168">' + "{0:.4f}".format(align.get_pi([pop])) + '</td>\n'
        html_string += '              </tr>\n'
html_string += '            </table>\n'
html_string += '          </td>\n'
html_string += '        </tr>\n'

# The between population differentiation section.
if pops != None and len(pops) > 1:
    html_string += '        <tr class="smallSpaceUnder">\n'
    html_string += '          <td style="font-size:30px; font-weight:bold">Between-population differentiation</td>\n'
    html_string += '        </tr>\n'
    html_string += '        <tr class="largeSpaceUnder">\n'
    html_string += '          <td style="font-family:helvetica; font-size:12px">\n'
    html_string += '            <table width="840" border="0" cellpadding="0" cellspacing="1">\n'
    html_string += '              <tr>\n'
    html_string += '                <td width="168" style="font-weight:bold">Population 1</td>\n'
    html_string += '                <td width="168" style="font-weight:bold">Population 2</td>\n'
    html_string += '                <td width="168" style="font-weight:bold"><a href="http://www.jstor.org/stable/2408641">F<sub>ST</sub></a></td>\n'
    html_string += '                <td width="168" style="font-weight:bold"><a href="http://onlinelibrary.wiley.com/doi/10.1111/mec.12842/abstract">d<sub>XY</a></td>\n'
    html_string += '                <td width="168" style="font-weight:bold">d<sub>f</td>\n'
    html_string += '              </tr>\n'
    for x in range(0,len(pops)-1):
        for y in range(x+1,len(pops)):
            html_string += '              <tr>\n'
            html_string += '                <td width="168" style="font-weight:bold">' + pops[x] + '</td>\n'
            html_string += '                <td width="168" style="font-weight:bold">' + pops[y] + '</td>\n'
            f_st = align.get_F_st([pops[x], pops[y]])
            if x == 0 and y == 1:
                html_string += '                <!-- First f_st: string start -->\n'
            if f_st is None:
                html_string += '                <td>NA</td>\n'
            else:
                html_string += '                <td width="168">' + "{0:.4f}".format(f_st) + '</td>\n'
            if x == 0 and y == 1:
                html_string += '                <!-- First f_st: string end -->\n'
            d_xy = align.get_d_xy([pops[x], pops[y]])
            if x == 0 and y == 1:
                html_string += '                <!-- First d_xy: string start -->\n'
            if d_xy is None:
                html_string += '                <td width="168">NA</td>\n'
            else:
                html_string += '                <td width="168">' + "{0:.4f}".format(d_xy) + '</td>\n'
            if x == 0 and y == 1:
                html_string += '                <!-- First d_xy: string end -->\n'
            d_f = align.get_d_f([pops[x], pops[y]])
            if x == 0 and y == 1:
                html_string += '                <!-- First d_f: string start -->\n'
            if d_f is None:
                html_string += '                <td width="168">NA</td>\n'
            else:
                html_string += '                <td width="168">' + "{0:.4f}".format(d_f) + '</td>\n'
            if x == 0 and y == 1:
                html_string += '                <!-- First d_f: string end -->\n'
            html_string += '              </tr>\n'
    html_string += '            </table>\n'
    html_string += '          </td>\n'
    html_string += '        </tr>\n'

# The genealogical sorting index section.
if len(pops) > 0:
    html_string += '        <tr class="smallSpaceUnder">\n'
    html_string += '          <td style="font-size:30px; font-weight:bold">Genealogical sorting index</td>\n'
    html_string += '        </tr>\n'
    html_string += '        <tr class="largeSpaceUnder">\n'
    html_string += '          <td style="font-family:helvetica; font-size:12px">\n'
    html_string += '            <table width="840" border="0" cellpadding="0" cellspacing="1">\n'
    html_string += '              <tr>\n'
    html_string += '                <td width="168" style="font-weight:bold">Population</td>\n'
    html_string += '                <td width="672" style="font-weight:bold"><a href="http://www.bioone.org/doi/abs/10.1111/j.1558-5646.2008.00442.x">gsi</a></td>\n'
    html_string += '              </tr>\n'
    for x in range(0,len(pops)):
        html_string += '              <tr>\n'
        html_string += '                <td width="168" style="font-weight:bold">' + pops[x] + '</td>\n'
        if x == 0:
            html_string += '                <!-- First gsi: string start -->\n'
        elif x == 1:
            html_string += '                <!-- Second gsi: string start -->\n'
        if gsis[x] == None:
            html_string += '                <td width="672">NA</td>\n'
        else:
            html_string += '                <td width="672">' + "{0:.4f}".format(gsis[x]) + '</td>\n'
        if x == 0:
            html_string += '                <!-- First gsi: string end -->\n'
        elif x == 1:
            html_string += '                <!-- Second gsi: string end -->\n'
        html_string += '              </tr>\n'
    html_string += '            </table>\n'
    html_string += '          </td>\n'
    html_string += '        </tr>\n'

# Finalize the html string.
html_string += '      </table>\n'
html_string += '    </div>\n'
html_string += '  </body>\n'
html_string += '</html>\n'
html_string += ''
html_string += ''

# Write the html string to STDOUT.
outfile.write(html_string)
