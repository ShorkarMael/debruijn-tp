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

__author__ = "Your Name"
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
            #print(line)
            yield line
    fastq.close()


def cut_kmer(sequence, k):
    for x in range(len(sequence) + 1 - k):
        kmer = sequence[x:x + k]
        yield kmer


def build_kmer_dict(fastq_file, k):
    kmer_dict = dict()
    for sequence in read_fastq(fastq_file):
        #print(sequence)
        for kmer in cut_kmer(sequence, k):
            kmer_dict[kmer] = kmer_dict.get(kmer, 0) + 1
    return kmer_dict


def build_graph(kmer_dict):
    graph = nx.DiGraph()
    for kmer, w in kmer_dict.items():
        prev_weight_int = 0
        k = len(kmer)
        #print(kmer[0:k-1], kmer[1:k], w)
        if graph.has_edge(kmer[0:k-1], kmer[1:k]):
            prev_weight = graph.get_edge_data(kmer[0:k-1], kmer[1:k], 'weight')
            prev_weight_int = prev_weight['weight']
        graph.add_edge(kmer[0:k-1], kmer[1:k], weight=w + prev_weight_int)
        #print(graph.edges(data=True))
        #draw_graph(graph, 'test.jpg')
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


def remove_paths(graph, path_list, delete_entry_node, delete_sink_node):
    pass

def std(data):
    pass


def select_best_path(graph, path_list, path_length, weight_avg_list, 
                     delete_entry_node=False, delete_sink_node=False):
    pass

def path_average_weight(graph, path):
    pass

def solve_bubble(graph, ancestor_node, descendant_node):
    pass

def simplify_bubbles(graph):
    pass

def solve_entry_tips(graph, starting_nodes):
    pass

def solve_out_tips(graph, ending_nodes):
    pass

def get_starting_nodes(graph):
    start_nodes = list()
    for node in graph.nodes:
        if nx.ancestors(graph, node) == set():
            start_nodes.append(node)
            #print('start : {}'.format(node))
    return start_nodes


def get_sink_nodes(graph):
    sink_nodes = list()
    for node in graph.nodes:
        if nx.descendants(graph, node) == set():
            sink_nodes.append(node)
            #print('sink : {}'.format(node))
    return sink_nodes

def get_contigs(graph, starting_nodes, ending_nodes):
    list_contigs = []
    for starting_node in starting_nodes:
        for ending_node in ending_nodes:
            for path in nx.all_simple_paths(graph, starting_node, ending_node):
                contig = [p[:-1] for p in path]
                contig.append(ending_node[-1])
                contig = ''.join(contig)
                list_contigs.append((contig,len(contig)))
    return list_contigs


def save_contigs(contigs_list, output_file):
    pass

#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    k = args.kmer_size
    tmp = 0
    graph = build_graph(build_kmer_dict(args.fastq_file, k))
    #print(graph.edges(data=True))
    #draw_graph(graph, 'test.jpg')
    get_starting_nodes(graph)
    get_sink_nodes(graph)

if __name__ == '__main__':
    main()
