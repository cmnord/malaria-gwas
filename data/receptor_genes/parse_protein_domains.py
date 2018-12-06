from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import sys
import os
import numpy as np
import pandas as pd

if len(sys.argv) < 2:
	print "you must call program as: python parse_alignments.py <protein_domains.txt> <sequences.fa>"
	sys.exit(1)

filename = sys.argv[1]
sequences = sys.argv[2]

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

def read_domains(domains):
	print "hi"
	name, dstats, p_doms = None, {}, {}
	cprim = None
	for line in domains:
		if not line.startswith('\t'):
			cprim = line[:-1]
			dstats[cprim] = {}
			p_doms[cprim] = {}
		elif line.startswith('\t') and not line.startswith('\t\t'):
			tm = line.find("numTM")
			rel = line.find("reliability")
			# print line[tm+7:tm+8]
			# print line[rel+13:-3]
			dstats[cprim]['numTm'] = line[tm+7:tm+8]
			dstats[cprim]["reliability"] = line[rel+13:-3]
		elif line.startswith('\t\t'):
			l = line.find("loc")
			loc = line[l+5:-4]
			f = line.find("from")
			fe = line[f+6:].find("\"")
			fr = line[f+6:f+6+fe]
			t = line.find("to")
			tt = line[t+4:].find("\"")
			to = line[t+4:t+4+tt]
			if loc == "O":
				if "O" in p_doms[cprim]:
					p_doms[cprim]["O"].append([fr, to])
				else:
					p_doms[cprim]["O"] = [[fr, to]]
			elif loc == "I":
				if "I" in p_doms[cprim]:
					p_doms[cprim]["I"].append([fr, to])
				else:
					p_doms[cprim]["I"] = [[fr, to]]
			else:
				if "M" in p_doms[cprim]:
					p_doms[cprim]["M"].append([fr, to])
				else:
					p_doms[cprim]["M"] = [[fr, to]]

	return dstats, p_doms


with open(filename) as fp:
	dstats, domains = read_domains(fp)
fp.close()




def make_outer_domain_fasta(domains, sequences):
	names = []
	reads = []
	with open(sequences) as f:
		for (name, seq) in read_fasta(f):
			names.append(name[1:])
			reads.append(seq)
	f.close()
	new_filename = sequences[:-3] + '_outer.fa'
	newfile = []
	for i in range(len(names)):
		header = '>' + names[i]
		newfile.append(header)
		try:
			outer_domains = domains[names[i]]["O"]
			seq = list(reads[i])[int(outer_domains[0][0]):int(outer_domains[0][1])]
			newfile.append("".join(seq))
		except: pass
	return new_filename, newfile

new_filename, newfile = make_outer_domain_fasta(domains, sequences)
with open(new_filename, 'w') as f:
	for line in newfile:
		f.write(line + '\n')




    #for (name, seq) in read_fasta(fp):
	 #    names.append(name[1:])
	 #    reads.append(seq)
	 #    nn = False
		# if name[1:] in infected:
		# 	infected_names.append(name[1:])
		# 	infected_reads.append(seq)
		# 	nn = True
	 #    if nn == False:
		# 	nf_names.append(name[1:])
		# 	nf_reads.append(seq)