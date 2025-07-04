
## read data to network
import networkx as nx
import random
import copy
import os
import pickle
import numpy
import csv

#!/usr/bin/env python
# -*- coding: utf-8 -*- 
# @Time    : 6.6.2019 14.53
# @Author  : YINYIN
# @Site    : 
# @File    : construct_network.py
# @Software: PyCharm

## read data to network
import networkx as nx
import random
import copy
import os
import pickle
import numpy
import csv

try:
    from scipy.stats import rankdata
except:
    print("scipy is not installed, rank-based distance methods won't work")

# Get the shortest path for those who have been calculated.
from functools import wraps

def memo_2(func):
    distance = {}
    @wraps(func)
    def _inner_wraps_2(G, source_id, target_id):
        if (source_id, target_id) in distance:
            result = distance[(source_id, target_id)]
        elif (target_id, source_id) in distance:
            result = distance[(target_id, source_id)]
        else:
            result = func(G, source_id, target_id)
            distance[(source_id, target_id)] = result
        return result
    return _inner_wraps_2

@memo_2
def get_shortest_path_length_between(G, source_id, target_id):
    return nx.shortest_path_length(G, source_id, target_id)

class Construct_Network:

    def __init__(self, file_name):
        self.filename = file_name
        self.setNode, self.setEdge, self.dictNode, self.dictEdge = self.get_nodes_and_edges_from_csv_file()
        self.g = self.create_network_from_sif_file(include_unconnected=True)
        self.G = self.get_network(only_lcc=True)

    def get_nodes_and_edges_from_csv_file(self):
        """
        Parse the CSV file into node and edge sets and dictionaries.
        """
        self.setNode = set()
        self.setEdge = set()  # Initialize to an empty set instead of None
        self.dictNode = {}
        self.dictEdge = {}
        with open(self.filename, 'r') as f:
            reader = csv.reader(f)
            next(reader)  # Skip header row
            for row in reader:
                if len(row) >= 2:
                    entrez_a = row[0]
                    entrez_b = row[1]
                    self.setNode.add(entrez_a)  # Add node A
                    self.setNode.add(entrez_b)  # Add node B
                    self.setEdge.add((entrez_a, entrez_b))  # Add edge between A and B
        
        if len(self.setEdge) == 0:  # Handle case where no edges are found
            print("Warning: No edges found in the file!")
            self.setEdge = set()  # Ensure it is an empty set, not None
        
        if len(self.dictNode) == 0: self.dictNode = None
        if len(self.dictEdge) == 0: self.dictEdge = None
        return self.setNode, self.setEdge, self.dictNode, self.dictEdge

    def create_network_from_sif_file(self, include_unconnected=True):
        """
        Use the node and edge data to create the network graph.
        """
        # Check if setNode is valid
        if not self.setNode:
            raise ValueError("self.setNode is None or empty. Please check the data loading process.")

        self.g = nx.Graph()  # Create an empty graph

        # If include_unconnected is True, add all nodes, even the unconnected ones
        if include_unconnected:
            self.g.add_nodes_from(self.setNode)

        # Check if setEdge is valid and add edges to the graph
        if self.setEdge:
            print(f"Adding edges: {self.setEdge}")  # Debug output
            self.g.add_edges_from(self.setEdge)
        else:
            print("Warning: self.setEdge is None or empty. No edges will be added to the graph.")
        
        # If no edges are added, consider adding logic to prevent an empty graph
        if len(self.g.nodes) == 0:
            raise ValueError("No nodes were added to the graph. Please check the input data.")
        
        # Return the constructed graph
        return self.g

    def get_network(self, only_lcc):
        """
        Returns the largest connected component (LCC) of the network if only_lcc is True
        """
        if only_lcc:
            print("Shrinking network to its LCC", len(self.g.nodes()), len(self.g.edges()))
            # Find connected components and return the largest one
            result_list = [c for c in sorted(nx.connected_components(self.g), key=len, reverse=True)]
            if not result_list:
                print("Warning: No connected components found.")
                self.G = self.g  # Return the original graph
            else:
                self.G = nx.subgraph(self.g, result_list[0])  # NetworkX subgraph method wrapper
                print("Final shape:", len(self.G.nodes()), len(self.G.edges()))
        else:
            self.G = self.g  # Return the original graph if no LCC filtering is needed
        return self.G

    def save_lcc(self):
        network_lcc_file = self.filename + ".lcc"
        if not os.path.exists(network_lcc_file):
            with open(network_lcc_file, 'w') as f:
                for u, v in self.G.edges():
                    f.write("%s 1 %s\n" % (u, v))

if __name__ == '__main__':
    file_name = "data/Interactome2022.csv"
    network = Construct_Network(file_name)
