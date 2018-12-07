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
p = sequences.find('\\')
protein = sequences[:p]
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

# read in the protein_domains file and find the sequences in heac domain
with open(filename) as fp:
	dstats, domains = read_domains(fp)
fp.close()

#print domains

# makes list that hass all domans, to be used for multiple alignment for each of the outer domains
def make_outer_domain_fasta_array(domains, sequences):
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
		#try:
		print names[i]
		outer_domains = domains[names[i]]["O"]
		o_d = []
		for j in range(len(outer_domains)):
			seq = list(reads[i])[int(outer_domains[j][0]):int(outer_domains[j][1])]
			o_d.append("".join(seq))
		newfile.append(o_d)
		#except: pass
	return new_filename, newfile

new_filename, newfile = make_outer_domain_fasta_array(domains, sequences)

# make new multiple alignment fastas for each outer domain
def make_ma_fasta(newfile, new_filename):
	for i in range(len(newfile[1])):
		with open(new_filename[:-3]+str(i+1)+new_filename[-3:] , 'w') as f:
			for line in newfile:
				if line[0] == '>':
					f.write(line + '\n')
				else:
					f.write(line[i] + '\n')
#make_ma_fasta(newfile, new_filename)

#make a folder with individual documents containing the human, and comparison sequence
def make_pairwise_fastas(protein, newfile):
	HUMAN_SEQ = None
	folder = "pairwise1"
	f1 = os.path.join(protein, folder)
	path = os.path.join(directory, f1)
	for i in range(len(newfile)):
		if newfile[i] == ">homo_sapiens":
			HUMAN_SEQ = newfile[i+1]
			break
	print HUMAN_SEQ

	for i in range(len(newfile)):
		if newfile[i] != ">homo_sapiens" and newfile[i] != "" and newfile[i][0] == ">":
			if newfile[i+1] == "":
				continue
			print newfile[i]
			fn = "hs_" + newfile[i][1:] + ".fa"
			with open(os.path.join(path,fn), 'w') as f:
				f.write(">homo_sapiens\n")
				f.write(HUMAN_SEQ + '\n')
				f.write(newfile[i]+'\n')
				f.write(newfile[i+1])

make_pairwise_fastas(protein, newfile)