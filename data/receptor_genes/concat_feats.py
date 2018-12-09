from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import sys
import os
import numpy as np
import pandas as pd


### Useful features
tfr1_IFs = pd.read_csv("tfr1/tfr1_impfeats_aanuc.csv")
scarb1_IFs = pd.read_csv("scarb1/outer/scarb1_impfeats_aanuc.csv")
duffy_aaIFs = pd.read_csv("duffy/outer/one_hot_important_aag.csv")
duffy_nucIFs = pd.read_csv("duffy/outer/duffy_impfeats_nucs.csv")

### One-hot encoded features
print "###############################################################################"



### This scrips pulls together all the useful features, and constructs a new feature matrix that uses only those.

tfr1_IFs = tfr1_IFs.drop(["Unnamed: 0","Overall"],axis=1)
scarb1_IFs = scarb1_IFs.drop(["Unnamed: 0","Overall"], axis=1)
#print scarb1_IFs
#print duffy_nucIFs
duffy_nucIFs = duffy_nucIFs.drop(["Unnamed: 0","Overall"],axis=1)
#print duffy_nucIFs
duffy_aaIFs = duffy_aaIFs.drop(["Unnamed: 0","Overall"],axis=1)
#print duffy_aaIFs

#print tfr1_IFs

allimps = pd.concat([tfr1_IFs, scarb1_IFs,duffy_aaIFs])
#print allimps

# use this to iterate over all of the useful features
for f in allimps['feats']:
	print f