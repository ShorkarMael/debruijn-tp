#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""Perform assembly based on debruijn graph."""

import argparse
import os
import sys
import networkx as nx
import matplotlib
from operator import itemgetter
import random

random.seed(9001)
from random import randint
import statistics
import matplotlib.pyplot as plt

__author__ = "Shorkar Maël"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Your Name"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Your Name"
__email__ = "your@email.fr"
__status__ = "Developpement"


def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
    "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', dest='fastq_file', type=isfile,
                        required=True, help="Fastq file")
    parser.add_argument('-k', dest='kmer_size', type=int,
                        default=21, help="K-mer size (default 21)")
    parser.add_argument('-o', dest='output_file', type=str,
                        default=os.curdir + os.sep + "contigs.fasta",
                        help="Output contigs in fasta file")
    return parser.parse_args()


def read_fastq(fastq_file):
    fastq = open(fastq_file, 'r')
    for lines in fastq:
        if lines[0] == 'A' or lines[0] == 'T' or lines[0] == 'C' or lines[0] == 'G':
            line = lines.strip('\n')
            # print(line)
            yield line
    fastq.close()


def cut_kmer(sequence, k):
    for x in range(len(sequence) + 1 - k):
        kmer = sequence[x:x + k]
        yield kmer


def build_kmer_dict(fastq_file, k):
    kmer_dict = dict()
    for sequence in read_fastq(fastq_file):
        # print(sequence)
        for kmer in cut_kmer(sequence, k):
            kmer_dict[kmer] = kmer_dict.get(kmer, 0) + 1
    return kmer_dict


def build_graph(kmer_dict):
    graph = nx.DiGraph()
    for kmer, w in kmer_dict.items():
        prev_weight_int = 0
        k = len(kmer)
        # print(kmer[0:k-1], kmer[1:k], w)
        if graph.has_edge(kmer[0:k - 1], kmer[1:k]):
            prev_weight = graph.get_edge_data(kmer[0:k - 1], kmer[1:k], 'weight')
            prev_weight_int = prev_weight['weight']
        graph.add_edge(kmer[0:k - 1], kmer[1:k], weight=w + prev_weight_int)
        # print(graph.edges(data=True))
        # draw_graph(graph, 'test.jpg')
    return graph


def draw_graph(graph, graphimg_file):
    """Draw the graph
    """
    fig, ax = plt.subplots()
    elarge = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] > 3]
    # print(elarge)
    esmall = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] <= 3]
    # print(elarge)
    # Draw the graph with networkx
    # pos=nx.spring_layout(graph)
    pos = nx.random_layout(graph)
    nx.draw_networkx_nodes(graph, pos, node_size=6)
    nx.draw_networkx_edges(graph, pos, edgelist=elarge, width=6)
    nx.draw_networkx_edges(graph, pos, edgelist=esmall, width=6, alpha=0.5,
                           edge_color='b', style='dashed')
    # nx.draw_networkx(graph, pos, node_size=10, with_labels=False)
    # save image
    plt.savefig(graphimg_file)


def std(data):
    return statistics.stdev(data)


def select_best_path(graph, path_list, path_length, weight_avg_list,
                     delete_entry_node=False, delete_sink_node=False):
    random.seed(9001)
    #  First get the best weight
    max_weight = max(weight_avg_list)
    best_weights = []
    #  List of best weight index
    for i in range(len(weight_avg_list)):
        if weight_avg_list[i] == max_weight:
            best_weights.append(i)
    if len(best_weights) > 1:
        max_length = max(path_length)
        best_of_best = []
        #  Best weight paths
        for i in best_weights:
            if path_length[i] >= max_length:
                max_length = path_length[i]
        for i in best_weights:
            if path_length[i] == max_length:
                best_of_best.append(i)
        if len(best_of_best) > 1:
            best_path = random.sample(best_of_best, 2) # random sample in population of the best
            best_path = best_path[0]
        else:
            best_path = best_of_best[0]
    else:
        best_path = best_weights[0]
    path_to_remove = path_list[:best_path] + path_list[best_path + 1:]
    graph = remove_paths(graph, path_to_remove, delete_entry_node, delete_sink_node)
    return graph


def path_average_weight(graph, path):
    list_weight = []
    for i, j, edge in graph.subgraph(path).edges(data=True):
        list_weight.append(edge['weight'])
    mean = statistics.mean(list_weight)
    return mean


def remove_paths(graph, path_list, delete_entry_node, delete_sink_node):
    for path in path_list:
        if delete_entry_node == True and delete_sink_node == False:
            graph.remove_nodes_from(path[:-1])
        elif delete_entry_node and delete_sink_node:
            graph.remove_nodes_from(path)
        elif delete_entry_node == False and delete_sink_node == True:
            graph.remove_nodes_from(path[1:])
        elif delete_entry_node == False and delete_sink_node == False:
            graph.remove_nodes_from(path[1:-1])
    return graph

def solve_bubble(graph, ancestor_node, descendant_node):
    # List of all simple paths
    list_paths = list(nx.all_simple_paths(graph, ancestor_node, descendant_node))
    path_length = []
    weight_avg_list = []
    for path in list_paths:
        path_length.append(len(path))
        weight_avg_list.append(path_average_weight(graph, path))
    graph = select_best_path(graph, list_paths, path_length, weight_avg_list)

    return graph


def simplify_bubbles(graph):
    for node in graph.nodes:
        predecessors = list(graph.predecessors(node))
        if len(predecessors) > 1:
            for i in range(len(predecessors)):
                for j in range(i, len(predecessors)):
                    path = nx.lowest_common_ancestor(graph, predecessors[i], predecessors[j]) #  getting common ancestors
                    if path != predecessors[i] and path != predecessors[j]:
                        graph = solve_bubble(graph, path, node)
                        return graph
    return graph # If nothing is simplified return the same graph


def solve_entry_tips(graph, starting_nodes):
    # First create a list of nodes with more than 1 ancestor
    target = []
    for node in graph.nodes():
        ancestor = list(graph.predecessors(node))
        if len(ancestor) > 1:
            target.append(node)
    # If no predecessors no modification
    if len(target) == 0:
        return graph
    list_paths = []
    weight_avg_list = []
    path_length = []
    for i in range(len(starting_nodes)):
        list_paths.append(list(nx.all_simple_paths(graph, starting_nodes[i], target[0])))
        list_paths[i] = list_paths[i][0]
        weight_avg_list.append(path_average_weight(graph, list_paths[i]))
        path_length.append(len(list_paths[i]))
    # Best path at the end
    graph = select_best_path(graph, list_paths, path_length, weight_avg_list, True, False)
    return graph


def solve_out_tips(graph, ending_nodes):
    # First create a list of nodes with more than 1 ancestor
    target = []
    for node in graph.nodes():
        ancestor = list(graph.descendants(node))
        if len(ancestor) > 1:
            target.append(node)
    # If no predecessors no modification
    if len(target) == 0:
        return graph
    list_paths = []
    weight_avg_list = []
    path_length = []
    for i in range(len(ending_nodes)):
        list_paths.append(list(nx.all_simple_paths(graph, ending_nodes[i], target[0])))
        list_paths[i] = list_paths[i][0]
        weight_avg_list.append(path_average_weight(graph, list_paths[i]))
        path_length.append(len(list_paths[i]))
    # Best path at the end
    graph = select_best_path(graph, list_paths, path_length, weight_avg_list, False, True)
    return graph


def get_starting_nodes(graph):
    start_nodes = list()
    for node in graph.nodes:
        if nx.ancestors(graph, node) == set():
            start_nodes.append(node)
            # print('start : {}'.format(node))
    return start_nodes


def get_sink_nodes(graph):
    sink_nodes = list()
    for node in graph.nodes:
        if nx.descendants(graph, node) == set():
            sink_nodes.append(node)
            # print('sink : {}'.format(node))
    return sink_nodes


def get_contigs(graph, starting_nodes, ending_nodes):
    list_contigs = []
    for starting_node in starting_nodes:
        for ending_node in ending_nodes:
            for path in nx.all_simple_paths(graph, starting_node, ending_node):
                contig = [p[:-1] for p in path]
                contig.append(ending_node[-1])
                contig = ''.join(contig)
                list_contigs.append((contig, len(contig)))
    return list_contigs


def save_contigs(contigs_list, output_file):
    with open(output_file, 'w+') as f:
        for i in range(len(contigs_list)):
            f.write('>contig_{} len={}\n'.format(i, contigs_list[i][1]))
            f.write(fill(contigs_list[i][0]))
            f.write('\n')
    f.close()


def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i + width] for i in range(0, len(text), width))


# ==============================================================
# Main program
# ==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    k = args.kmer_size
    tmp = 0
    graph = build_graph(build_kmer_dict(args.fastq_file, k))
    # print(graph.edges(data=True))
    # draw_graph(graph, 'test.jpg')
    start = get_starting_nodes(graph)
    end = get_sink_nodes(graph)
    list_contig = get_contigs(graph,start,end)
    save_contigs(list_contig,'test_save.fasta')
    removed_graph = remove_paths(graph, ['TCAGAGCTCTAGAGTTGGTT','TCAATCACACCCACCACGTG'], True, False).edges(data=True)


if __name__ == '__main__':
    main()
