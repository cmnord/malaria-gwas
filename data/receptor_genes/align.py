import sys
import os

if len(sys.argv) < 1:
    print("you must call program as: python align.py <directory with seqs>")
    sys.exit(1)

directory = sys.argv[1]
print(directory)
#file2 = sys.argv[2]

dirr = os.getcwd()
print(dirr)

workingdir = os.path.join(dirr, directory)


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


reads2 = []
names2 = []


for filename in os.listdir(workingdir):
    reads = []
    names = []

    if "homo_sapiens" in filename:
        u1 = filename.find('_')+1
        prot = filename[:u1-1]
        spec = filename[u1:-3]
        #u2 = spec.find('_')
        #spec = spec[:u2]

        with open(os.path.join(workingdir, filename)) as fp:
            for (name, seq) in read_fasta(fp):
                names.append(name)
                reads.append(seq)
        fp.close()

        with open(os.path.join(workingdir, prot + '_all_primates.fa'), 'a') as f3:
            f3.write('>' + spec + '\n')
            for read in reads:
                f3.write(read)
        f3.close()
        continue

for filename in os.listdir(workingdir):
    reads = []
    names = []

    if filename.endswith(".fa"):
        u1 = filename.find('_')+1
        prot = filename[:u1-1]
        spec = filename[u1:-3]
        if "homo_sapiens" in spec or "all_primates" in spec:
            continue
        with open(os.path.join(workingdir, filename)) as fp:
            for (name, seq) in read_fasta(fp):
                names.append(name)
                reads.append(seq)
        fp.close()
        with open(os.path.join(workingdir, prot + '_all_primates.fa'), 'a') as f3:
            f3.write('\n' + '>' + spec + '\n')
            for i in range(len(reads)):
                f3.write(reads[i])
        f3.close()


# ### GLOBAL ALIGNMENTS -- all best
# #alignments = pairwise2.align.globalxx(X,Y)

# ### LOCAL ALIGNMENTS -- all best
# #alignments = pairwise2.align.localxx(X,Y)

# ### GLOBAL WITH GAP PENALTY -- all best
# ### 2 match, -1 mismatch, -0.5 open gap, -0.1 extend gap
# alignments = pairwise2.align.globalms(X,Y,5,-4,-2,-0.5)
# print(format(alignments[0], 'psl'))

# # for a in alignments:
# # 	print((format_alignment(*a)))

# print(len(alignments))


# from Bio.Align.Applications import ClustalwCommandline
# cline = ClustalwCommandline("clustalw2", infile="opuntia.fasta")
# import os
