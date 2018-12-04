from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import sys
import os

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

#X = "ACGGGT"
#Y = "ACG"

# filename = sys.argv[1]
# filename2 = sys.argv[2]
# else:
#     sys.exit("Usage: " + sys.argv[0] + " [fasta_filename]")
reads2 = []
names2 = []


directory = 'C:/Users/Kari/Documents/MIT/Senior/Fall 2018/6.047/malaria-gwas/data/receptor_genes/tfr1/protein'

for filename in os.listdir(directory):
	reads = []
	names = []

	if "homo_sapiens" in filename:
		u1 = filename.find('_')+1
		prot = filename[:u1-1]
		spec = filename[u1:-3]
		#u2 = spec.find('_')
		#spec = spec[:u2]

		with open(filename) as fp:
		    for (name, seq) in read_fasta(fp):
		        names.append(name)
		        reads.append(seq)
		fp.close()

		with open(prot + '_all_primates.fa', 'a') as f3:
			f3.write('>' + spec + '\n')
			for read in reads:
				f3.write(read)
		continue

for filename in os.listdir(directory):
	reads   = []
	names   = []

	if filename.endswith(".fa") and filename.startswith("tfr1"):
		u1 = filename.find('_')+1
		prot = filename[:u1-1]
		spec = filename[u1:-3]
		#u2 = spec.find('_')
		#spec = spec[:u2]
		if "homo_sapiens" in spec or "all_primates" in spec:
			continue
		with open(filename) as fp:
		    for (name, seq) in read_fasta(fp):
		        names.append(name)
		        reads.append(seq)
		fp.close()
		# with open(prot + '_test2_p1_names.fa', 'a') as f:
		# 	f.write('>' + spec + '\n')# + ' ' + names[0][1:])
		# 	f.write(reads[0] + '\n')
		# f.close()
		# with open(prot + '_test2_p2_names.fa', 'a') as f2:
		# 	if len(reads) >= 2:
		# 		f2.write('>' + spec + '\n')# + ' ' + names[1][1:])
		# 		f2.write(reads[1] + '\n')
		# f2.close()
		with open(prot + '_all_primates.fa', 'a') as f3:
			f3.write('\n' + '>' + spec + '\n')
			for i in range(len(reads)):
				f3.write(reads[i])


# with open('band3_allseq_p1.fa') as fp:
#     for(name, seq) in read_fasta(fp):
#         names2.append(name)
#         reads2.append(seq)
# fp.close()


# X = reads[0]
# Y = reads2[0]

# #print reads
# #print names2

# ### GLOBAL ALIGNMENTS -- all best
# #alignments = pairwise2.align.globalxx(X,Y)

# ### LOCAL ALIGNMENTS -- all best
# #alignments = pairwise2.align.localxx(X,Y)

# ### GLOBAL WITH GAP PENALTY -- all best
# ### 2 match, -1 mismatch, -0.5 open gap, -0.1 extend gap
# alignments = pairwise2.align.globalms(X,Y,5,-4,-2,-0.5)
# print format(alignments[0], 'psl')

# # for a in alignments:
# # 	print (format_alignment(*a))

# print len(alignments)


# from Bio.Align.Applications import ClustalwCommandline
# cline = ClustalwCommandline("clustalw2", infile="opuntia.fasta")
# import os
