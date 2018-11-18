import pandas as pd 
import numpy as np 
import csv
import statistics
import matplotlib.pyplot as plt
import argparse
import os

def get_parser():
	parser = argparse.ArgumentParser()
	parser.add_argument('--f',
		type=str,
		help='Screening data with gene locations and strand. Include full filepath')
	parser.add_argument('--add1',
		type=int,
		default=4,
		help='Number of nucleotides to examine past the 5 prime end of the seq')
	parser.add_argument('--add2',
		type=int,
		default=6,
		help='Number of nucleotides to examine past the 3 prime end of the seq')
	return parser


if __name__ == '__main__':
	args = get_parser().parse_args()
	print "FILENAME: " + args.f
	alldata = pd.read_csv(args.f)

	e = args.f.find('.csv')
	name = args.f[:e]
	if args.add1 == 0:
		name = name+'_proto'
	print 'BED FILE NAME IS: '+name+'.bed'

	with open(name+'.bed', 'w') as fi:
		for guide in alldata.index:
			try:
				seq = alldata.loc[guide]['Location']
				strand = alldata.loc[guide]['Strand']
				gid = alldata.loc[guide]['ID']

				colon = seq.find(':')
				chrom = seq[0:colon]
				hyph = seq.find('-')

				start = seq[colon+1:hyph]
				end = seq[hyph+1:]			

				if strand == '+':
					tmerstart = str(int(start)-args.add1)
					tmerend = str(int(end)+args.add2)
				else:
					tmerstart = str(int(start)-args.add2)
					tmerend = str(int(end)+args.add1)

				full = chrom + '\t' + tmerstart + '\t' + tmerend + '\t' + gid + '\t' + '1' + '\t' + strand + '\n'

				fi.write(full)
			except: 
				raise Exception("The datafile is missing some needed data. make sure \'Location\', \'Strand\', and \'ID\' present for all gudies")