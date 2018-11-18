import pandas as pd 
import numpy as np 
import csv
import matplotlib.pyplot as plt
import argparse
import os

def get_parser():
	parser = argparse.ArgumentParser()
	parser.add_argument('--f',
		type=str,
		help='Screening data with gene locations and strand. Include full filepath')
	return parser


if __name__ == '__main__':
	args = get_parser().parse_args()
	print "FILENAME: " + args.f
	alldata = pd.read_csv(args.f)

	e = args.f.find('.csv')
	name = args.f[:e]
	print 'BED FILE NAME IS: '+name+'.bed'

	with open(name+'.bed', 'w') as fi:
		for guide in alldata.index:
			seq = int(alldata.loc[guide]['Position'])
			end = seq+1
			SNP_name = alldata.loc[guide]['SNP']
			chrom = alldata.loc[guide]['Chromosome']

			full = '{}\t{}\t{}\t{}\n'.format(chrom, seq, end, SNP_name)
			fi.write(full)
