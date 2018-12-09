from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import sys
import os
import numpy as np
import pandas as pd


### Useful features
# tfr1_IFs = pd.read_csv("tfr1/tfr1_impfeats_aanuc.csv")
# scarb1_IFs = pd.read_csv("scarb1/outer/scarb1_impfeats_aanuc.csv")
# duffy_aaIFs = pd.read_csv("duffy/outer/one_hot_important_aag.csv")
# duffy_nucIFs = pd.read_csv("duffy/outer/duffy_impfeats_nucs.csv")

# ### One-hot encoded features
# print "###############################################################################"



# ### This scrips pulls together all the useful features, and constructs a new feature matrix that uses only those.

# tfr1_IFs = tfr1_IFs.drop(["Unnamed: 0","Overall"],axis=1)
# scarb1_IFs = scarb1_IFs.drop(["Unnamed: 0","Overall"], axis=1)
# #print scarb1_IFs
# #print duffy_nucIFs
# duffy_nucIFs = duffy_nucIFs.drop(["Unnamed: 0","Overall"],axis=1)
# #print duffy_nucIFs
# duffy_aaIFs = duffy_aaIFs.drop(["Unnamed: 0","Overall"],axis=1)
# #print duffy_aaIFs

# #print tfr1_IFs

# allimps = pd.concat([tfr1_IFs, scarb1_IFs,duffy_aaIFs])
# #print allimps

# # use this to iterate over all of the useful features

# for f in allimps['feats']:
# 	print f


folders = ["gypc","basigin","scarb1","tfr1","epha2"]
numout = [1, 1, 1, 1, 1]
directory = os.getcwd()
filedict = {}
outernum = 1
for i in range(len(folders)):
	folder = folders[i]
	for onehot in ["","aagroups"]:
		outernum = 1
		while outernum <= numout[i] and outernum!=None:
			print folder, outernum, numout[i]
			doc = pd.read_csv(folder+"\outer\\"+folder+"_all_primates_outer"+str(outernum)+onehot+".csv")
			if len(doc.columns.values) > 2:
				filedict[folder+"_"+str(outernum)+"_"+onehot] = doc
				outernum+=1
			else:
				outernum=None
infected = ['pan_troglodytes', 'gorilla_gorilla', 'pongo_abelii', 'pongo_pygmaeus', 'homo_sapiens', 'aotus_nancymaae', 'aotus trivirgatus', 'saimiri_boliviensis']


species = ["saguinus_labiatus", "pongo_pygmaeus" ,             
"nomascus_leucogenys","pan_troglodytes",           
"pan_paniscus","homo_sapiens",     
"gorilla_gorilla","rhinopithecus_roxellana",  
"piliocolobus_tephrosceles","macaca_mulatta",              
"chlorocebus_sabaeus","cercocebus_atys",        
"mandrillus_leucophaeus","papio_anubis",                
"galeopterus_variegatus","plecturocebus_moloch",        
"aotus_nancymaae","aotus_trivirgatus",           
"callithrix_jacchus","colobus_angolensis ",
"macaca_fascicularis","rhinopithecus_bieti",         
"theropithecus_gelada","macaca_nemestrina",           
"chlorocebus_aethiops","pongo_abelii"]

d = {'Unnamed: 0': species}

totfeats = pd.DataFrame(data=d)
for key in filedict:
	spec = filedict[key]
	spec = spec.drop(["Species"], axis=1)
	totfeats = pd.merge(totfeats,spec, how='inner',on='Unnamed: 0')

print totfeats

inf = []
for spec in totfeats["Unnamed: 0"]:
	if spec in infected:
		inf.append("infected")
	else: inf.append("resistant")
totfeats["Species"] = inf

print totfeats

df1 = totfeats.loc[:,(totfeats!=0).any()]
df2 = df1.loc[:,(df1!=1).any()]
df3 = df2.loc[:, (df2.sum(axis=0))>1]
print df3
df4 = df3.loc[:, (df3.sum(axis=0))<6]
df4["Species"] = df3["Species"]
df4["Unnamed: 0"] = totfeats["Unnamed: 0"]
print df4

df4.to_csv("all_combined_features_minband3duff.csv")


# for key in filedict:
# 	print len(filedict[key].columns.values)
# 	print filedict[key]
# 	print key