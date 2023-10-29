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
import textwrap
import matplotlib.pyplot as plt
matplotlib.use("Agg")

__author__ = "Lebib Ines"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Lebib Ines"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Lebib Ines"
__email__ = "leines74@gmail.com"
__status__ = "Developpement"

def isfile(path): # pragma: no cover
    """Check if path is an existing file.

    :param path: (str) Path to the file
    
    :raises ArgumentTypeError: If file doesn't exist
    
    :return: (str) Path 
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments(): # pragma: no cover
    """Retrieves the arguments of the program.

    :return: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', dest='fastq_file', type=isfile,
                        required=True, help="Fastq file")
    parser.add_argument('-k', dest='kmer_size', type=int,
                        default=22, help="k-mer size (default 22)")
    parser.add_argument('-o', dest='output_file', type=str,
                        default=os.curdir + os.sep + "contigs.fasta",
                        help="Output contigs in fasta file (default contigs.fasta)")
    parser.add_argument('-f', dest='graphimg_file', type=str,
                        help="Save graph as an image (png)")
    return parser.parse_args()


def read_fastq(fastq_file):
    """Extract reads from fastq files.
    :param fastq_file: (str) Path to the fastq file.
    :return: A generator object that iterate the read sequences. 
    """
    with open(fastq_file, 'r') as fastq :
        for line in fastq :
            yield(next(fastq).strip())
            next(fastq)
            next(fastq)



def cut_kmer(read, kmer_size):
    read_length = len(read)
    if read_length < kmer_size:
        print("le k_mer est plus petit que la longueur du read indiqué")
    
    for i in range(read_length - kmer_size + 1):
        kmer = read[i:i + kmer_size]
        yield kmer


def build_kmer_dict(fastq_file, kmer_size) :
    """Build a dictionnary object of all kmer occurrences in the fastq file

    :param fastq_file: (str) Path to the fastq file.
    :return: A dictionnary object that identify all kmer occurrences.
    """
    kmer_dict = {}
    for seq in read_fastq(fastq_file) : 
        for kmer in cut_kmer(seq, kmer_size) : 
            if kmer in kmer_dict : 
                kmer_dict[kmer] += 1 
            else : 
                kmer_dict[kmer] = 1
    
    return kmer_dict


def build_graph(kmer_dict: dict[str, int]) :
    """Build the debruijn graph

    :param kmer_dict: A dictionnary object that identify all kmer occurrences.
    :return: A directed graph (nx) of all kmer substring and weight (occurrence).
    """
    graph = nx.DiGraph()

    # A chaque iteration, la clé kmer et l'occurence sont assignes
    # aux variables
    
    for kmer, occurrence in kmer_dict.items():
        # Recupere  tous les caracteres du début jusqu'au deuxieme
        # caracteres avec la fin
        prefix = kmer[:-1]

        # Recupere tous les caracteres à partir du deuxieme 
        # caractere jusqu'à la fin
        suffix = kmer[1:]
        
        if prefix not in graph:
            graph.add_node(prefix)
        if suffix not in graph:
            graph.add_node(suffix)
        
        graph.add_edge(prefix, suffix, weight=occurrence)
    
    return graph


def remove_paths(graph, path_list, delete_entry_node, delete_sink_node):
    """Remove a list of path in a graph. A path is set of connected node in
    the graph

    :param graph: (nx.DiGraph) A directed graph object
    :param path_list: (list) A list of path
    :param delete_entry_node: (boolean) True->We remove the first node of a path 
    :param delete_sink_node: (boolean) True->We remove the last node of a path
    :return: (nx.DiGraph) A directed graph object
    """
    remove_graph = graph.copy()

    for path in path_list:
        if delete_sink_node & delete_entry_node:
             for node in path:
                remove_graph.remove_node(node)
             break
        if delete_entry_node:
            first_node = path[0]
            remove_graph.remove_node(first_node)

        if delete_sink_node:
            last_node = path[-1]
            remove_graph.remove_node(last_node)
        else:
            for node in path[1:-1]:
                remove_graph.remove_node(node)
            

    return remove_graph

def select_best_path(graph, path_list, path_length, weight_avg_list, 
                     delete_entry_node=False, delete_sink_node=False):
    """Select the best path between different paths

    :param graph: (nx.DiGraph) A directed graph object
    :param path_list: (list) A list of path
    :param path_length_list: (list) A list of length of each path
    :param weight_avg_list: (list) A list of average weight of each path
    :param delete_entry_node: (boolean) True->We remove the first node of a path 
    :param delete_sink_node: (boolean) True->We remove the last node of a path
    :return: (nx.DiGraph) A directed graph object
    """
    stdev_w = statistics.stdev(weight_avg_list)
   
    stdev_l = statistics.stdev(path_length)
   
    path_ind = None
   
    if stdev_w > 0:
        path_ind = weight_avg_list.index(max(weight_avg_list))
    elif stdev_l > 0:
        path_ind = path_length.index(max(path_length))
    else:
        path_ind = random.randint(0, len(path_length) - 1)

    path_list.pop(path_ind)
    modified_graph = graph.copy()

    modified_graph = remove_paths(graph, path_list , delete_entry_node, delete_sink_node)
   
    return modified_graph


def path_average_weight(graph, path):
    """Compute the weight of a path

    :param graph: (nx.DiGraph) A directed graph object
    :param path: (list) A path consist of a list of nodes
    :return: (float) The average weight of a path
    """
    return statistics.mean([d["weight"] for (u, v, d) in graph.subgraph(path).edges(data=True)])

def solve_bubble(graph, ancestor_node, descendant_node):
    """Explore and solve bubble issue

    :param graph: (nx.DiGraph) A directed graph object
    :param ancestor_node: (str) An upstream node in the graph 
    :param descendant_node: (str) A downstream node in the graph
    :return: (nx.DiGraph) A directed graph object
    """
    paths = list(nx.all_simple_paths(graph, ancestor_node, descendant_node))
    weigths = [path_average_weight(graph, path) for path in paths]
    lengths = [len(path) - 1 for path in paths]
    solved_bubble = select_best_path(graph, paths, lengths, weigths)
    return solved_bubble


def simplify_bubbles(graph):
    """Detect and explode bubbles

    :param graph: (nx.DiGraph) A directed graph object
    :return: (nx.DiGraph) A directed graph object
    """
    bubble = False
    for node in graph.nodes:
        predecessors = list(graph.predecessors(node))

        if len(predecessors) > 1:
            for i in range(len(predecessors)):
                for j in range(i + 1, len(predecessors)):
                    ancestor_node = nx.lowest_common_ancestor(graph, predecessors[i], predecessors[j])

                    if ancestor_node is not None:
                        bubble = True
                        break

        if bubble:
            break

    if bubble:
        ancestor = ancestor_node
        node = node
        graph = solve_bubble(graph, ancestor, node)
        return simplify_bubbles(graph)

    return graph


def solve_entry_tips(graph, starting_nodes):
    """Remove entry tips

    :param graph: (nx.DiGraph) A directed graph object
    :return: (nx.DiGraph) A directed graph object
    """
    for node in graph:
        node_pred = list(graph.predecessors(node))
        if len(node_pred) > 1:
            valid_paths = []
            for node_start_i in starting_nodes:
                paths = list(nx.all_simple_paths(graph, node_start_i, node))
                valid_paths.extend([path[0] for path in paths if len(path) > 0])
            
            lengths = [len(path) - 1 for path in valid_paths]
            weights = [path_average_weight(graph, path) if lengths[i] > 1 else graph.get_edge_data(*path)["weight"]
                    for i, path in enumerate(valid_paths)]

            graph = select_best_path(graph, valid_paths, lengths, weights, delete_entry_node=True, delete_sink_node=False)
            graph = solve_entry_tips(graph, starting_nodes)
            break

    return graph

def solve_out_tips(graph, ending_nodes):
    """Remove out tips

    :param graph: (nx.DiGraph) A directed graph object
    :return: (nx.DiGraph) A directed graph object
    """
    for node in graph:
        node_success = list(graph.successors(node))
        if len(node_success) > 1:
            valid_paths = []
            for node_end_i in ending_nodes:
                paths = list(nx.all_simple_paths(graph, node, node_end_i))
                valid_paths.extend([path[0] for path in paths if len(path) > 0])
            
            lengths = [len(path) - 1 for path in valid_paths]
            weights = [path_average_weight(graph, path) if lengths[i] > 1 else graph.get_edge_data(*path)["weight"]
                    for i, path in enumerate(valid_paths)]

            graph = select_best_path(graph, valid_paths, lengths, weights, delete_entry_node=False, delete_sink_node=True)
            graph = solve_out_tips(graph, ending_nodes)
            break

    return graph

def get_starting_nodes(graph):
    """Get nodes without predecessors

    :param graph: (nx.DiGraph) A directed graph object
    :return: (list) A list of all nodes without predecessors
    """
    nodes_without_predecessors = [node for node in graph.nodes() if not any(graph.predecessors(node))]
    return nodes_without_predecessors

def get_sink_nodes(graph):
    """Get nodes without successors

    :param graph: (nx.DiGraph) A directed graph object
    :return: (list) A list of all nodes without successors
    """
    nodes_without_successors = []
    for node in graph.nodes():
        if not any(graph.successors(node)):
            nodes_without_successors.append(node)
    return nodes_without_successors


def get_contigs(graph, starting_nodes, ending_nodes):
    """Extract the contigs from the graph

    :param graph: (nx.DiGraph) A directed graph object 
    :param starting_nodes: (list) A list of nodes without predecessors
    :param ending_nodes: (list) A list of nodes without successors
    :return: (list) List of [contiguous sequence and their length]
    """
    contigs = []

    for start_node in starting_nodes:
        for end_node in ending_nodes:
            if nx.has_path(graph, start_node, end_node):
                chemins = list(nx.all_simple_paths(graph, start_node, end_node))
                contig = chemins[0][0]
                for path in chemins[0][1:]:
                    contig += path[-1]
            longueur_contig = len(contig)
            contigs.append((contig, longueur_contig))

    return contigs

def save_contigs(contigs_list, output_file):
    """Write all contigs in fasta format

    :param contig_list: (list) List of [contiguous sequence and their length]
    :param output_file: (str) Path to the output file
    """
    with open(output_file, 'w') as file:
        for i, (contig, length) in enumerate(contigs_list):
            file.write(f'>contig_{i} len={length}\n')
            enveloped_contig = textwrap.fill(contig, width=80)
            file.write(enveloped_contig)
            file.write('\n')


def draw_graph(graph, graphimg_file): # pragma: no cover
    """Draw the graph

    :param graph: (nx.DiGraph) A directed graph object
    :param graphimg_file: (str) Path to the output file
    """                                   
    fig, ax = plt.subplots()
    elarge = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] > 3]
    #print(elarge)
    esmall = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] <= 3]
    #print(elarge)
    # Draw the graph with networkx
    #pos=nx.spring_layout(graph)
    pos = nx.random_layout(graph)
    nx.draw_networkx_nodes(graph, pos, node_size=6)
    nx.draw_networkx_edges(graph, pos, edgelist=elarge, width=6)
    nx.draw_networkx_edges(graph, pos, edgelist=esmall, width=6, alpha=0.5, 
                           edge_color='b', style='dashed')
    #nx.draw_networkx(graph, pos, node_size=10, with_labels=False)
    # save image
    plt.savefig(graphimg_file)


#==============================================================
# Main program
#==============================================================
def main(): # pragma: no cover
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()

    fastq_file = args.fastq_file  # Mettez le chemin vers votre fichier FASTQ ici
    kmer_size = args.kmer_size
    kmer_dict = build_kmer_dict(fastq_file, kmer_size)
    debruijn_graph = build_graph(kmer_dict)
    starting_nodes = get_starting_nodes(debruijn_graph)
    sink_nodes = get_sink_nodes(debruijn_graph)
    contigs = get_contigs(debruijn_graph, starting_nodes, sink_nodes)
    output_file = args.output_file
    save_contigs(contigs, output_file)
    # Voir le dictionnaire construit
    # Fonctions de dessin du graphe
    # A decommenter si vous souhaitez visualiser un petit
    # graphe
    # Plot the graph
    if args.graphimg_file:
        draw_graph(debruijn_graph, args.graphimg_file)




if __name__ == '__main__': # pragma: no cover
    main()
