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

infected = ['pan_troglodytes', 'gorilla', 'pongo_abelii', 'homo_sapiens']
sometimes_infected = ['aotus_nancymaae', 'aotus trivirgatus', 'saimiri_boliviensis']

names = []
reads = []
infected_reads = []
infected_names = []
nf_reads = []
nf_names = []
prim_order = []



with open(filename) as fp:
    for (name, seq) in read_fasta(fp):
        names.append(name[1:])
        reads.append(seq)
        nn = False
        for n in infected:
        	if n in name:
        		infected_names.append(name[1:])
        		infected_reads.append(seq)
        		prim_order.append(n)
        		nn = True
        if nn == False:
    		nf_reads.append(seq)
    		nf_names.append(name[1:])

fp.close()

reads2 = []
for i in range(len(reads)):
	reads2.append(list(reads[i]))
aligndf = pd.DataFrame(reads2, columns=range(1,len(reads[0])+1), index=names)
aligndf.to_csv('tfr1/primates/tfr1_test3_alignment.csv')
print aligndf

# #print reads[0]

# nuc_dict = [None for i in range(len(reads[0]))]
# inf_dict = [None for i in range(len(reads[0]))]
# nf_dict = [None for i in range(len(reads[0]))]
# #print nuc_dict
# nuc_count = []
# inf_count = []
# nf_count = []
# nuc_max = []
# inf_max = []
# nf_max = []

# for i in range(len(nuc_dict)):
# 	nuc_dict[i] = {}
# 	for seq in reads:
# 		if seq[i] != '-':
# 			if seq[i] not in nuc_dict[i]:
# 				nuc_dict[i][seq[i]] = 1
# 			else:
# 				nuc_dict[i][seq[i]] += 1
# 	maxNuc = None
# 	maxCount = 0
# 	for key in nuc_dict[i]:
# 		if nuc_dict[i][key] > maxCount:
# 			maxNuc = key
# 			maxCount = nuc_dict[i][key]
# 	nuc_count.append(maxCount)
# 	nuc_max.append(maxNuc)

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

# assert len(nuc_count) == len(inf_count)

# compare = np.column_stack((nuc_count, inf_count, nf_count, nuc_max, inf_max, nf_max))
# print compare

# headers = [i for i in range(len(nuc_count))]

# compdf = pd.DataFrame(compare, columns=['all', 'inf', 'not_inf', 'all_nuc', 'inf_nuc', 'not_inf_nuc'])
# print compdf

# compdf.to_csv('tfr1\primates\compare1.csv')
