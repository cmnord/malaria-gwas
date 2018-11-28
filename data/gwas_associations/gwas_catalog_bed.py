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
	print alldata

	dup_dict = {}
	with open(name+'.bed', 'w') as fi:
		for guide in alldata.index:
			seq = alldata.loc[guide]['CHR_POS']
			SNP_name = alldata.loc[guide]['MAPPED_GENE']
			chrom = alldata.loc[guide]['CHR_ID']
			if chrom != 'X' and chrom != 'Y':
				try:
					chrom = int(chrom)
					end = int(seq)+1
				except Exception as e:
					print(e)
					continue
			if chrom == 'X':
				end = int(seq)+1

			if chrom not in dup_dict:
				dup_dict[chrom] = [seq]
			else:
				if seq not in dup_dict[chrom]:
					dup_dict[chrom].append(seq)
				else:
					continue
			full = 'chr{}\t{}\t{}\t{}\n'.format(chrom, seq, end, SNP_name)
			fi.write(full)

		
			
