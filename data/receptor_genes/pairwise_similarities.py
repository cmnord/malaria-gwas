from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import sys
import os
import numpy as np
import pandas as pd

##################################
### TAKES IN THE OUTPUT OF PAIRWISE.SH AND ANNOTATES
###  AMINO ACIDS THAT ARE IDENTICAL TO THE HUMAN AMINO ACID
### WRITE A TXT TILE THAT ALLOWS FOR COMPARISON

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

infected = ['pan_troglodytes', 'gorilla_gorilla', 'pongo_abelii', 'pongo_pygmaeus', 'homo_sapiens', 'aotus_nancymaae', 'aotus trivirgatus', 'saimiri_boliviensis']

names = []
reads = []
infected_reads = []
infected_names = []
nf_reads = []
nf_names = []

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

def find_identical_subsequences(reads, use_blosum = False):
	ident = []
	for i in range(len(reads[0])):
		if reads[1][i] == '-' or reads[0][i] == '-':
			ident.append(0)
		elif reads[1][i] in blosum2[reads[0][i]]:
			ident.append(1)
		else: ident.append(0)
	return ident

if __name__ == "__main__":
	if len(sys.argv) < 1:
		print "you must call program as: python pairwise_similarities.py <alignment.afa>"
		sys.exit(1)

	filename = sys.argv[1]

	directory = os.getcwd()

	with open(filename) as fp:
	    for (name, seq) in read_fasta(fp):
	        names.append(name[1:])
	        reads.append(seq)
	fp.close()

	HUMAN = None
	ident_dict = {}
	# find human sequences
	for i in range(len(names)):
		if names[i] == "homo_sapiens":
			HUMAN = reads[i]
	for i in range(len(names)):
		if names[i] != "homo_sapiens":
			compare = [HUMAN,reads[i]]
			n = ["homo_sapiens",names[i]]
			ident_dict[n[1]] = find_identical_subsequences(compare, True)

	for key in ident_dict:
		print key, ident_dict[key]


fn = os.path.join(directory, filename[:-4])
with open(fn+'.txt', 'w') as f:
	for key in ident_dict:
		f.write(key + '\t')
		for el in list(ident_dict[key]):
			f.write(str(el)+'\t')
		f.write('\n')