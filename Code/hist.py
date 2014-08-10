import pandas as pd
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import itertools
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.spatial.distance import seuclidean 
import scipy
# Import all needed libraries.

df = pd.read_csv('../Data/dataset.csv', header = 0)
AA_percents = {}
# Create a dictionary{} that will hold the amino acid names as a key and 
#the percent of each amino acid letter in the protein's AA sequence
Distances = []

for row in df.index:
    key = df.ix[row, 'viral_protein_name'] + ' , ' + df.ix[row, 'virus_species']
    value = ProteinAnalysis(df.ix[row, 'viral_protein_AA_seq']).get_amino_acids_percent()
    AA_percents[key] = pd.DataFrame.from_dict(value, orient='index')
# Find all possible pairs of proteins
# With scipy's euclidean distance function find the euclidean distannce between the AA seq percent
for x,y in itertools.combinations(AA_percents.values(), 2):
    Distances.append(scipy.spatial.distance.euclidean(x,y))
# Create a histogram that plots the euclidean distances of the proteins
numBins = 125
plt.hist(Distances, numBins,color='green')
plt.title("Amino Acid Comparison")
plt.xlabel("Distances between AA_seqs")
plt.ylabel("Frequencey ")
plt.savefig('Amino Acid Comparison')






    
