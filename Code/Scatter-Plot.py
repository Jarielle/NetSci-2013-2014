import pandas as pd
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import itertools
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.spatial.distance import seuclidean 
import scipy
import networkx as nx
from pylab import *
# Import all needed libraries.

df = pd.read_csv('../Data/dataset.csv', header = 0)
AA_percents = {}
# Create a dictionary{} that will hold the amino acid names as a key and 
#the percent of each amino acid letter in the protein's AA sequence
Distances = []
for row in df.index:
    key = df.ix[row, 'viral_protein_name'] + ' : ' + df.ix[row, 'virus_species']
    value = ProteinAnalysis(df.ix[row, 'viral_protein_AA_seq']).get_amino_acids_percent()
    AA_percents[key] = pd.DataFrame.from_dict(value, orient='index')
Distances = {}
# Find all possible pairs of proteins
# With scipy's euclidean distance function find the euclidean distannce between the AA seq percent
for x,y in itertools.combinations(AA_percents.keys(), 2):
    Distances[x,y] = scipy.spatial.distance.euclidean(AA_percents[x], AA_percents[y])

my_list = []
pos_edges = []
e069 = []

for x,y in Distances.items():
    for z in x:
        my_list.append(z)
    PE = x + (y,)
    pos_edges.append(PE)
# Create a list that holds all possible edges within the network. 
proteins = list(set(my_list))
# Create a list that holds all proteins.
for edge in pos_edges:
    for weight in edge:
        if .000<= weight <.069:
           e069.append(edge)
        else:
            pass
G = nx.DiGraph()
# Creates a Undirected graph
G.add_nodes_from(proteins)
# Adds the list proteins[] as nodes in the network.
pos=nx.spring_layout(G)
G.add_weighted_edges_from(e069, weight='weight')
#  Adds the list e069[] as edges in the network.
G_ud = G.to_undirected()
# Find the closeness and betweenness centrality for each node in the network.
x = nx.betweenness_centrality(G).values()
y = nx.closeness_centrality(G).values()
fig = plt.figure(figsize=(5, 5))
# Create a scatter plot that plots both the closeness and betweenness centrality for each node in the network.
ax = fig.add_subplot(1,1,1) 
ax.set_title("Amini Acid Comparison")
ax.set_xlabel("Betweenness Centrality")
ax.set_ylabel("Closeness Centrality")
