from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import sys
import os
import numpy as np
import pandas as pd

if len(sys.argv) < 1:
	print "you must call program as: python parse_alignments.py <alignment.afa>"
	sys.exit(1)

filename = sys.argv[1]

directory = os.getcwd()

def read_fasta(fp):
    name, seq = None, []
    for line in fp:
        line = line.rstrip()
        if line.startswith(">"):
            if name: yield (name, ''.join(seq))
            name, seq = line, []
        else:
            seq.append(line)
    if name: yield (name, ''.join(seq))

#filename = 'C:/Users/Kari/Documents/MIT/Senior/Fall 2018/6.047/malaria-gwas/data/receptor_genes/tfr1/primates/tfr1_test3_full.afa'

infected = ['pan_troglodytes', 'gorilla', 'pongo_abelii', 'homo_sapiens', 'aotus_nancymaae', 'aotus trivirgatus', 'saimiri_boliviensis']

names = []
reads = []
infected_reads = []
infected_names = []
nf_reads = []
nf_names = []

# read in the sequence, store the names and the corresponding reads
# nothing done with this piece of code right now
# maybe want to find sequences where all nucleotides are identical
# - these
with open(filename) as fp:
    for (name, seq) in read_fasta(fp):
        names.append(name[1:])
        reads.append(seq)
        nn = False
    	if name[1:] in infected:
    		infected_names.append(name[1:])
    		infected_reads.append(seq)
    		nn = True
        if nn == False:
    		nf_names.append(name[1:])
    		nf_reads.append(seq)

fp.close()

inf_matrix = []
for i in range(len(infected_reads)):
	inf_matrix.append(list(infected_reads[i]))

inf_r_np = np.array(inf_matrix)
inf_r_np_t = inf_r_np.transpose()

unique_inf = []
for i in range(len(inf_r_np_t)):
	unique_inf.append(np.unique(inf_r_np_t[i]))
unique_inf = np.array(unique_inf)
#unique_inf = np.unique(inf_r_np_t)
ident = []
print inf_r_np
print inf_r_np_t
print len(inf_r_np_t)
#print unique_inf

# find sequences in infected species that are identical
for i in range(len(unique_inf)):
	if len(unique_inf[i]) == 1 and unique_inf[i] != '-':
		ident.append(1)
	else:
		ident.append(0)

print ident
print len(ident)


# make a dataframe and csv file with the alignments
def make_alignment_csv():
	reads2 = []
	for i in range(len(reads)):
		reads2.append(list(reads[i]))
	aligndf = pd.DataFrame(reads2, columns=range(1,len(reads[0])+1), index=names)
	aligndf.to_csv(directory + '\\' + filename[:-4] + '_alignment.csv')
	print aligndf

# aa_dict = [None for i in range(len(reads[0]))]
# inf_dict = [None for i in range(len(reads[0]))]
# nf_dict = [None for i in range(len(reads[0]))]

# aa_count = []	# had the number of species containing the most common nucleotide
# inf_count = []	# has the number of infected containing the msot common nucleotide
# nf_count = []
# aa_max = []
# inf_max = []
# nf_max = []

# # make a dictionary for each position and add all amino acids found at that position
# for i in range(len(aa_dict)):
# 	aa_dict[i] = {}
# 	for seq in reads:
# 		if seq[i] != '-':
# 			if seq[i] not in aa_dict[i]:
# 				aa_dict[i][seq[i]] = 1
# 			else:
# 				aa_dict[i][seq[i]] += 1
# 	maxNuc = None
# 	maxCount = 0
# 	for key in aa_dict[i]:
# 		if aa_dict[i][key] > maxCount:
# 			maxNuc = key
# 			maxCount = aa_dict[i][key]
# 	aa_count.append(maxCount)
# 	aa_max.append(maxNuc)

# #which_prim = []

# for i in range(len(inf_dict)):
# 	inf_dict[i] = {}
# 	for j in range(len(infected_reads)):
# 		if infected_reads[j][i] != '-':
# 			if infected_reads[j][i] not in inf_dict[i]:
# 				inf_dict[i][infected_reads[j][i]] = 1
# 			else:
# 				inf_dict[i][infected_reads[j][i]] += 1
# 	maxNuc = None
# 	maxCount = 0
# 	for key in inf_dict[i]:
# 		if inf_dict[i][key] > maxCount:
# 			maxNuc = key
# 			maxCount = inf_dict[i][key]

# 	inf_count.append(maxCount)
# 	inf_max.append(maxNuc)

# for i in range(len(nf_dict)):
# 	nf_dict[i] = {}
# 	for j in range(len(nf_reads)):
# 		if nf_reads[j][i] != '-':
# 			if nf_reads[j][i] not in nf_dict[i]:
# 				nf_dict[i][nf_reads[j][i]] = 1
# 			else:
# 				nf_dict[i][nf_reads[j][i]] += 1
# 	maxNuc = None
# 	maxCount = 0
# 	for key in nf_dict[i]:
# 		if nf_dict[i][key] > maxCount:
# 			maxNuc = key
# 			maxCount = nf_dict[i][key]

# 	nf_count.append(maxCount)
# 	nf_max.append(maxNuc)

# assert len(aa_count) == len(inf_count)

# compare = np.column_stack((aa_count, inf_count, nf_count, aa_max, inf_max, nf_max))
# print compare

# headers = [i for i in range(len(aa_count))]

# compdf = pd.DataFrame(compare, columns=['all', 'inf', 'not_inf', 'all_nuc', 'inf_nuc', 'not_inf_nuc'])
# print compdf

# compdf.to_csv('tfr1\primates\compare1.csv')
