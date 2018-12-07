from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import sys
import os
import numpy as np
import pandas as pd

blosum = {"A":['A', 'C', 'G', 'S', 'T', 'V'],
			"R":['N', 'R', 'Q', 'E', 'H', 'K'],
			"N":['R', 'N', 'D', 'Q', 'E', 'G', 'H', 'K', 'S', 'T'],
			"D":['N', 'D', 'Q', 'E', 'S'],
			"C":['A', 'C'],
			"Q":['R', 'N', 'D', 'Q', 'E', 'H', 'K', 'M', 'S'],
			"E":['R', 'N', 'D', 'Q', 'E', 'H', 'K', 'S'],
			"G":['A', 'N', 'G', 'S'],
			"H":['R', 'N', 'Q', 'E', 'H', 'Y'],
			"I":['I', 'L', 'M', 'F', 'V'],
			"L":['I', 'L', 'M', 'F', 'V'],
			"K":['R', 'N', 'Q', 'E', 'K', 'S'],
			"M":['Q', 'I', 'L', 'M', 'F', 'V'],
			"F":['I', 'L', 'M', 'F', 'W', 'Y'],
			"P":['P'],
			"S":['A', 'N', 'D', 'Q', 'E', 'G', 'K', 'S', 'T'],
			"T":['A','N', 'S', 'T', 'V'],
			"W":['F', 'W', 'Y'],
			"Y":['H', 'F', 'W', 'Y'],
			"V":['A', 'I', 'L', 'M', 'T', 'V']}

blosum2 = {"A":['A','S'],
			"R":['R', 'Q', 'K'],
			"N":['N', 'D', 'S'],
			"D":['N', 'D','E'],
			"C":['C'],
			"Q":['R', 'Q', 'E', 'K'],
			"E":['D', 'Q', 'E','K'],
			"G":['G'],
			"H":['N', 'H', 'Y'],
			"I":['I', 'L', 'M', 'V'],
			"L":['I', 'L', 'M', 'V'],
			"K":['R', 'Q', 'E', 'K'],
			"M":['I', 'L', 'M', 'V'],
			"F":['F', 'W', 'Y'],
			"P":['P'],
			"S":['A', 'N', 'S', 'T'],
			"T":['S', 'T'],
			"W":['F', 'W', 'Y'],
			"Y":['H', 'F', 'W', 'Y'],
			"V":['I', 'L', 'M', 'V']}

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

infected = ['pan_troglodytes', 'gorilla_gorilla', 'pongo_abelii', 'homo_sapiens', 'aotus_nancymaae', 'aotus trivirgatus', 'saimiri_boliviensis']

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

def find_identical_subsequences(reads, use_blosum = False):
	read_matrix = []
	for i in range(len(reads)):
		read_matrix.append(list(reads[i]))

	reads_np = np.array(read_matrix)
	reads_np_t = reads_np.transpose()

	unique = []
	for i in range(len(reads_np_t)):
		unique.append(np.unique(reads_np_t[i]))
	unique = np.array(unique)
	ident = []
	# find sequences in infected species that are identical
	for i in range(len(unique)):
		if len(unique[i]) == 1:
			ident.append(1)
		else:
			if use_blosum == True:
				cons = True
				for j in range(len(reads_np_t[i])):
					for k in range(len(reads_np_t[i])):
						if reads_np_t[i][k] == '-' or reads_np_t[i][j] == '-':
							cons = False
						elif reads_np_t[i][k] not in blosum2[reads_np_t[i][j]]:
							cons = False
				if cons == True:
					ident.append(1)
				else:
					ident.append(0)
			else:
				ident.append(0)

	return ident

# make a dataframe and csv file with the alignments
# endname is something like '_alignments.csv' or whatever you want the file to end with
# m_file is true of false, depending on if you want a csv file
def make_alignment_csv(reads, endname, m_file):
	reads2 = []
	for i in range(len(reads)):
		reads2.append(list(reads[i]))
	if m_file == True:
		aligndf = pd.DataFrame(reads2, columns=range(1,len(reads[0])+1))
		aligndf.to_csv(os.path.join(directory,filename[:-4]) + endname)
	return reads2


all_ident = find_identical_subsequences(reads)
inf_ident = find_identical_subsequences(infected_reads)
res_ident = find_identical_subsequences(nf_reads)
ident_list = [all_ident, inf_ident, res_ident]

all_identB = find_identical_subsequences(reads, True)
inf_identB = find_identical_subsequences(infected_reads, True)
res_identB = find_identical_subsequences(nf_reads, True)
ident_listB = [all_identB, inf_identB, res_identB]

# print ident_list
ident_matrix = np.array(make_alignment_csv(ident_list, '_b2_identical.csv', False))


# find all sequences of 1's
def find_conserved_seqs(seqs):
	pairs = []
	i = 1
	prev = seqs[0]
	if prev == 1:
		tup = [0, None]
	else:
		tup = [None, None]
	while i < len(seqs):
		if seqs[i] == 1 and prev == 1:
			prev = 1
			i+=1
		elif seqs[i] == 1 and prev == 0 :
			tup = [i, None]
			prev = 1
			i+=1
		elif seqs[i] == 0 and prev == 1:
			tup[1] = i-1
			if tup[1]-tup[0] <= 5:
				pass
			else:
				pairs.append(tup)
			tup = [None, None]
			i+=1
			prev = 0
		else:
			tup = [None, None]
			prev = 0
			i+=1
	return pairs

inf_pairs = find_conserved_seqs(inf_ident)
res_pairs = find_conserved_seqs(res_ident)
all_pairs = find_conserved_seqs(all_ident)
print inf_pairs
print res_pairs
print all_pairs

inf_pairsB = find_conserved_seqs(inf_identB)
res_pairsB = find_conserved_seqs(res_identB)
all_pairsB = find_conserved_seqs(all_identB)


# find parts that are similarly conserved between two species
# donn'''ttttt think this acutally says anything significant.
# will probably need to compare the conserved parts among infected species
#	with each other species independently. for each animal that doens't 
#	get infected by vivax, we could align it with the consensus alignment,
#	or the sequences that are completely conserved in the infected primates
#	
def find_overlaps(pairs1, pairs2):
	overlaps = []
	if len(pairs1) < len(pairs2):
		first = pairs2
		second = pairs1
	else:
		first = pairs1
		second = pairs2
	for i in range(len(first)):
		pair1 = first[i]
		for j in range(len(second)):
			pair2 = second[j]
			if pair2[0] in range(pair1[0]-2, pair1[1]) and pair2[1] in range(pair1[0], pair1[1]+2):
				overlaps.append([pair1, pair2])
			elif pair1[0] in range(pair2[0]-2, pair2[1]) and pair1[1] in range(pair2[0], pair2[1]+2):
				overlaps.append([pair1, pair2])
			elif pair2[1] < pair1[0]:
				pass
			elif pair2[0] > pair1[1]:
				continue
	return overlaps

inf_all_overlaps = find_overlaps(inf_pairs, all_pairs)
inf_all_overlapsB = find_overlaps(inf_pairsB, all_pairsB)
print np.array(inf_all_overlaps)
print np.array(inf_all_overlapsB)

# for l in inf_all_overlaps:
# 	print infected_reads[0][l[0][0]:l[0][1]]
# 	print 

# for i in range(len(infected_reads[0])):
# 	if 



	# inf_dict = [None for i in range(len(reads[0]))]
	# nf_dict = [None for i in range(len(reads[0]))]
	# inf_count = []	# has the number of infected containing the msot common nucleotide
	# nf_count = []
	# inf_max = []
	# nf_max = []



def count_identical(specs):
	# takes in a list of lists of alignments
	# calculates the number of aa's in each position that are the same
	# finds the most common nucleotide for each position

	specs_dict = [None for i in range(len(reads[0]))]

	aa_count = []	# has the number of species containing the most common nucleotide
	aa_max = []		# has the nucleotide that is most frequently seen

	# make a dictionary for each position and add all amino acids found at that position
	for i in range(len(specs_dict)):
		specs_dict[i] = {}
		for seq in specs:
			if seq[i] != '-':
				if seq[i] not in specs_dict[i]:
					specs_dict[i][seq[i]] = 1
				else:
					specs_dict[i][seq[i]] += 1
		maxNuc = None
		maxCount = 0
		for key in specs_dict[i]:
			if specs_dict[i][key] > maxCount:
				maxNuc = key
				maxCount = specs_dict[i][key]
		aa_count.append(maxCount)
		aa_max.append(maxNuc)

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
