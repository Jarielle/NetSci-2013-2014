import pandas as pd
import matplotlib.pyplot as plt
#Import all needed python libraries.
ongc_list = []
#Create a list that will hold oncogenicity viruses
beng_list = []
#Create a list that will hold benign viruses.
df = pd.read_csv('dataset.csv', header = 0)
for row in df.index:
    if df.ix[row, 'virus_species_oncogenic'] == 0:
       beng_list.append(len(df.ix[row, 'viral_protein_AA_seq']))
    elif df.ix[row, 'virus_species_oncogenic'] == 1:
         ongc_list.append(len(df.ix[row, 'viral_protein_AA_seq']))
#Using pandas to set conditions that will added the protein to the ongc[] or beng[] list.  
data = [ongc_list,beng_list]
plt.boxplot(data)
#Plot the data variable in a box plot.
#Create a box plot that compares the length of oncogenic viruses and benign viruses.
labels = ('Oncogenic', 'Non-Oncogenic')
plt.legend()
plt.xticks(range(1,3),labels, rotation=15)
plt.xlabel('Oncogenicity')
plt.ylabel('Length of the Amino Acid Sequence')
plt.title('The Length of the amino acid sequence in oncogenic and non-oncogenic viral proteins')
plt.show()
