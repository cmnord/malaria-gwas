from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import sys
import os
import numpy as np
import pandas as pd

infected = ['pan_troglodytes', 'gorilla_gorilla', 'pongo_abelii', 'pongo_pygmaeus', 'homo_sapiens', 'aotus_nancymaae', 'aotus trivirgatus', 'saimiri_boliviensis']

map_aa = {"A":1,
			"R":2,
			"N":3,
			"D":4,
			"C":5,
			"Q":6,
			"E":7,
			"G":8,
			"H":9,
			"I":10,
			"L":11,
			"K":12,
			"M":13,
			"F":14,
			"P":15,
			"S":16,
			"T":17,
			"W":18,
			"Y":19,
			"V":20,
			'-':21}
aas = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V", "-"]
#aas = ["Q", "E", "H", "-"]

# chemsim = {"A":["A","G"], "R":["R","K"], "N":["N","D"], "D":["N","D"], "C":["C","S"], "Q":["N","Q","E","H"],
#             "E":["D","Q","E"], "G":["A","G"], "H":["Y","F","Q","H"], "I":["V","I","M","L"],
#             "L":["I","L","M","F"], "K":["R","K"], "M":["I","L","M","F"], "F":["H","L","F","M","W","Y"],
#             "P":["P"], "S":["C","S","T"], "T":["S","T"], "W":["F","W","Y"], "Y":["H","F","W","Y"], "V":["I","V"]}
chemsim = {0:["A","G"],
			1:["R","K"],
			2:["N","D","Q","E"],
			3:["C","S","T"],
			4:["I","L","V","M"],
			5:["P"],
			6:["Y","W","F","H"],
			7:["-"]}


print len(aas)
names = []
reads = []

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

def encode_aa(aa):
	aa_code = [0 for i in range(21)]
	index = map_aa[aa]-1
	aa_code[index] = 1
	return aa_code

def encode_groups(aa):
	aa_code = [0 for i in range(8)]
	for key in chemsim:
		if aa in chemsim[key]:
			aa_code[key] = 1
	return aa_code


if __name__ == "__main__":
	if len(sys.argv) < 1:
		print "you must call program as: python one_hot.py <alignment.afa>"
		sys.exit(1)

	filename = sys.argv[1]
	p = filename.find('\\')
	o = filename.find('.afa')
	protein = filename[:p]
	region = filename[o-6:o]

	directory = os.getcwd()

	with open(filename) as fp:
	    for (name, seq) in read_fasta(fp):
	        names.append(name[1:])
	        reads.append(seq)
	fp.close()

	# figure out how to name the features
	length = len(reads[0])*21
	#headers = [protein+"_"+region+"_"+str(int(i/21))+aas[i%21] for i in range(len(reads[0]*21))]
	headers = [protein+"_"+region+"_pos"+str(int(i/8))+'_'+str(i%8) for i in range(len(reads[0]*8))]

	#print headers[:30]
	# for k in range(len(headers)):
	# 	headers[k] = headers[k]+"_"+str(int(i/21))+aas[length%21]

	features = [[] for i in range(len(names))]
	for i in range(len(names)):
		p_feat = []
		for aa in reads[i]:
			#p_feat+=encode_aa(aa)
			p_feat+=encode_groups(aa)
		features[i]+= p_feat

	specs = []
	for name in names:
		if name in infected:
			specs.append("infected")
		else: specs.append("resistant")
	df = pd.DataFrame(features, columns=headers, index=names)
	df["Species"] = specs

	print df
	df.to_csv(os.path.join(directory,filename[:o])+'aagroups.csv')

