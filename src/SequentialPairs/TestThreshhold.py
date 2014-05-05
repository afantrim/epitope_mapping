from __future__ import division

'''
.. module:: TestThreshhold
   :platform: Unix, Windows
   :synopsis: This will be used to automatically generate Precision-Recall plots for varying MIN_ST values.

.. moduleauthor:: Amelia F. Antrim <amelia.f.antrim@gmail.com>

'''
from PDBUtil import PDBUtil
from PrecisionRecall import PrecisionRecall
import operator
from Statistics import Statistics
from Bio.PDB import *
import sys
import numpy
import networkx as nx

class TestThreshhold:

    @staticmethod
    def filter_components(est_pos, Pairs):
        graph = TestThreshhold.make_graph(est_pos, Pairs)
        graph = TestThreshhold.get_components(graph)
        for pos in est_pos:
            if not graph.has_node(pos):
                est_pos.remove(pos)
        return est_pos

    @staticmethod
    def get_components(graph, cutoff=2):
        # Note: chose 3 because functional epitopes are usually 3-5 residues in size
        components = nx.connected_component_subgraphs(graph)
        for component in components:
            if component.number_of_nodes() < cutoff:
                # Remove all the nodes in this subgraph from the graph overall
                nodes = component.nodes()
                for node in nodes:
                    graph.remove_node(node)
        return graph

    @staticmethod
    def make_graph(est_pos, Pairs):
        # Initialize a graph
        G = nx.Graph()
        # Consider the estimated positives to be vertices
        # Note: residues are named by their ID number and nothing else; most convert back to residue objects
        for residue in est_pos:
            G.add_node(residue)
        for pair in Pairs:
            if pair[0] in est_pos and pair[1] in est_pos:
                G.add_edge(pair[0], pair[1])
        # Their appearances in pairs will be the edges
        return G
    
    @staticmethod
    def evalutate_true_positives(co_structure, antigen_chain):
         # Get all the chains
        chains = []
        for chain in co_structure:
            if chain == antigen_chain:
                continue
            else:
                chains.append(chain)
        
        # Compare the antigen chain to all the other chains
        for chain1 in chains:
            if chain1 == antigen_chain:
                continue
            else:
                # Calculate the distance matrix
                # Note that the first element is the antigen chain so the chain1 dictionary
                # will be nested inside of that
                dist_tree = PDBUtil.calc_dist_matrix(antigen_chain, chain1)
                
                # Get a list of all residues in the chain
                antigen_list = Selection.unfold_entities(antigen_chain, 'R')
                i = 0
                for antigen_res in dist_tree.keys():
                    for antibody_res in dist_tree[antigen_res].keys():
                        if dist_tree[antigen_res][antibody_res] < dist and antigen_res not in true_positives:
                            i += 1
                            # The true positives is a list of antigen residues that are part of the epitope
                            true_positives.append(antigen_res)
        return true_positives
        

    @staticmethod
    def make_plot(contact_dict, co_structure, antigen_chain, Pairs, increment=1.0, dist=8.0):
        # Sort the contact dictionary to a list of tuples so that it can be iterated through more quickly
        est_pos = []
        precision = []
        recall = []
        
        true_positives = evaluate_true_positives(co_structure, antigen_chain)
    
        # Initialize the precision-recall plot
        plot = PrecisionRecall()
        
        print "Evaluting precision and recall..."
        
        threshhold = 1.0
        # Iterate through all the prescribed threshhold values
        while threshhold:
        
            # Reset the positives
            est_pos = []
            
            # Now get the estimated positives from the contact dictionary given this value
            for residue, count in contact_dict.iteritems():
                if count > threshhold:
                    est_pos.append(residue)
            
            # Make a graph to evaluate the interconnectedness of the residues
            est_pos = TestThreshhold.filter_components(est_pos, Pairs)
            
            # If there are no true or false positives at this threshhold, quit and make the plot
            if len(est_pos) == 0:
                break
            
            # Now iterate through the connected components of this graph and take out the ones that are too small
            
            stats = Statistics(est_pos, true_positives, antigen_chain)
            TP = stats.true_positives
            FP = stats.false_positives
            
            print "Threshhold:"
            print threshhold
            print "True Positives:"
            print true_positives
            print "Estimated True Positives:"
            print TP
    
            # Precision = true est pos/ true pos
            this_precision = (float(len(TP)) / float(len(true_positives)))
            precision.append(this_precision)
            
            # Recall = true est pos/(true pos + false pos)
            this_recall = (float(len(TP))/(float(len(TP)) + float(len(FP))))
            recall.append(this_recall)
            
            # Add this threshhold to the plot
            #plot.add_threshhold(precision, recall)
            
            # Increment the threshhold by the indicated amount
            threshhold += increment
        
        print "The maximum threshhold is %.2f." %(threshhold-increment)
        print "Done."
        
        # Make the plot
        plot.make_plot(precision, recall)
            