from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import sys
import os
import numpy as np
import pandas as pd

if len(sys.argv) < 1:
    print "you must call program as: python parse_alignments.py <sequences.fa>"
    sys.exit(1)

sequences = sys.argv[1]
p = sequences.find('outer')
protein = sequences[:p+5]
filename = sequences[p+6:]
directory = os.getcwd()
print protein, filename


def read_fasta(fp):
    name, seq = None, []
    for line in fp:
        line = line.rstrip()
        if line.startswith(">"):
            if name:
                yield (name, ''.join(seq))
            name, seq = line, []
        else:
            seq.append(line)
    if name:
        yield (name, ''.join(seq))


def make_pairwise_fastas(protein, names, reads):
    HUMAN_SEQ = None
    folder = "pairwise1"
    f1 = os.path.join(protein, folder)
    path = os.path.join(directory, f1)
    for i in range(len(reads)):
        if names[i] == ">homo_sapiens":
            HUMAN_SEQ = reads[i]
            break
    print HUMAN_SEQ

    for i in range(len(reads)):
        if names[i] != ">homo_sapiens" and reads[i] != "":
            if reads[i] == "":
                continue
            fn = "hs_" + names[i][1:] + ".fa"
            with open(os.path.join(path, fn), 'w') as f:
                f.write(">homo_sapiens\n")
                f.write(HUMAN_SEQ + '\n')
                f.write(names[i]+'\n')
                f.write(reads[i])


names = []
reads = []
with open(os.path.join(directory, sequences)) as fp:
    for (name, seq) in read_fasta(fp):
        names.append(name)
        reads.append(seq)

make_pairwise_fastas(protein, names, reads)
