from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import matplotlib.pyplot as plt
import numpy as np

fasta_file = 'sample.fasta'
subset_file = 'subset.txt'

'''takes a Uniprot proteome and a file containing a subset of Uniprot IDs from 
the same proteome and then plots the fraction of cysteines, as well as acid, 
basic and hydrophobic amino acids in the whole proteome and the subset.'''
    
proteome_cysteins = 0
proteome_basic_aa = 0
proteome_acidic_aa = 0
proteome_hydrophobic_aa = 0
subset_cysteins = 0
subset_basic_aa = 0
subset_acidic_aa = 0
subset_hydrophobic_aa = 0
    
full_sequence_length = 0
subset_sequence_length = 0
    
handle = open(fasta_file, 'r')
    
# calculating the number cysteins, basic amino acids, acidic amino acids, 
# and hydrophobic amino acids for the whole proteome
for record in SeqIO.parse(handle, 'fasta'):
    X = ProteinAnalysis(str(record.seq))
        
    print(record.id)
        
    proteome_cysteins += X.count_amino_acids()['C']
    # basic amino acids
    proteome_basic_aa += X.count_amino_acids()['R']
    proteome_basic_aa += X.count_amino_acids()['H']
    proteome_basic_aa += X.count_amino_acids()['K']
    # acidic amino acids
    proteome_acidic_aa += X.count_amino_acids()['D']
    proteome_acidic_aa += X.count_amino_acids()['E']
    # hydrophobic amino acids
    proteome_hydrophobic_aa += X.count_amino_acids()['A']
    proteome_hydrophobic_aa += X.count_amino_acids()['I']
    proteome_hydrophobic_aa += X.count_amino_acids()['L']
    proteome_hydrophobic_aa += X.count_amino_acids()['M']
    proteome_hydrophobic_aa += X.count_amino_acids()['F']
    proteome_hydrophobic_aa += X.count_amino_acids()['V']
    proteome_hydrophobic_aa += X.count_amino_acids()['P']
    proteome_hydrophobic_aa += X.count_amino_acids()['G']
    # append the length of sequence for later calculations
    full_sequence_length += len(record.seq)
        
handle.close()
        
print('Whole proteome: ', proteome_cysteins, 'cysteins,', proteome_basic_aa, 'basic aminos acids,', proteome_acidic_aa, 'acidic amino acids, and', proteome_hydrophobic_aa, 'hydrphobic amino acids')
print('The full length of the sequences is ', full_sequence_length)

# finding percent compositions of different types of amino acids in the whole proteome
fraction_cysteins = proteome_cysteins/full_sequence_length
fraction_basic = proteome_basic_aa/full_sequence_length
fraction_acidic = proteome_acidic_aa/full_sequence_length
fraction_hydrophobic = proteome_hydrophobic_aa/full_sequence_length
    
# opening subset file and extracting ids from each line, then adding them
# to a list
subset = open(subset_file, 'r')
subset_IDs = []
for line in subset:
    stripped_line = line.strip()
    subset_IDs.append(stripped_line)
subset.close()
    
# calculating the number cysteins, basic amino acids, acidic amino acids, 
# and hydrophobic amino acids for the subset
handle = open(fasta_file, 'r')
for record in SeqIO.parse(handle, 'fasta'):
    for seq_id in subset_IDs:
        if record.id == str(seq_id):
            X = ProteinAnalysis(str(record.seq))
            subset_cysteins += X.count_amino_acids()['C']
            # basic amino acids
            subset_basic_aa += X.count_amino_acids()['R']
            subset_basic_aa += X.count_amino_acids()['H']
            subset_basic_aa += X.count_amino_acids()['K']
            # acidic amino acids
            subset_acidic_aa += X.count_amino_acids()['D']
            subset_acidic_aa += X.count_amino_acids()['E']
            # hydrophobic amino acids
            subset_hydrophobic_aa += X.count_amino_acids()['A']
            subset_hydrophobic_aa += X.count_amino_acids()['I']
            subset_hydrophobic_aa += X.count_amino_acids()['L']
            subset_hydrophobic_aa += X.count_amino_acids()['M']
            subset_hydrophobic_aa += X.count_amino_acids()['F']
            subset_hydrophobic_aa += X.count_amino_acids()['V']
            subset_hydrophobic_aa += X.count_amino_acids()['P']
            subset_hydrophobic_aa += X.count_amino_acids()['G']
            # append the length of sequence for later calculations
            subset_sequence_length += len(record.seq)
handle.close()
                
print('Subset: ', subset_cysteins, 'cysteins,', subset_basic_aa, 'basic amino acids,', subset_acidic_aa, 'acidic amino acids, and', subset_hydrophobic_aa, 'hydrophobic amino acids')
print('The full length of the subset sequences is', subset_sequence_length)
    
# finding percent compositions of different types of amino acids in the subset
fraction_cysteins_subset = subset_cysteins/subset_sequence_length
fraction_basic_subset = subset_basic_aa/subset_sequence_length
fraction_acidic_subset = subset_acidic_aa/subset_sequence_length
fraction_hydrophobic_subset = subset_hydrophobic_aa/subset_sequence_length
    
# creating axis values
objects = ('Cystein (proteome)', 'Basic (proteome)', 'Acidic (proteome)', 'Hydrophobic (proteome)', 'Cystein (subset)', 'Basic (subset)', 'Acidic (subset)', 'Hydrophobic (subset)')
y_pos = np.arange(len(objects))
percentages = [fraction_cysteins,fraction_basic,fraction_acidic,fraction_hydrophobic,fraction_cysteins_subset,fraction_basic_subset,fraction_acidic_subset,fraction_hydrophobic_subset]
    
# plot graph
plt.bar(y_pos, percentages, align='center', alpha=0.5)
plt.xticks(y_pos, objects, rotation=90)
plt.ylabel('Percentages')
plt.xlabel('Amino Acids Type')
plt.title('Amino Acid Percentages in Whole Proteome versus the Subset')
    
plt.show()
