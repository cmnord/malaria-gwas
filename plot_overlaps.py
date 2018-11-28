import sys
import pandas as pd
import matplotlib.patches as mpatches
from matplotlib import pyplot as plt

def get_num(elem):
        for i in range(len(elem)):
            try:
                int(elem[i])
                return int(elem[i:])
            except:
                continue 

def filter_cell_types(meta):
    exluded_groups = {'iPSC'}
    exluded_types = {'fetal', 'culture', 'carcinoma', 'leukemia'}

    filtered_cell_types = []
    types_to_check = meta['Epigenome ID (EID)'].values[2:]
    groups_to_check = meta['GROUP'].values[2:]
    descripts_to_check = meta['Standardized Epigenome name'].values[2:]

    for i in range(len(types_to_check)):
        valid = True
        for word in exluded_groups: # filters by group
            if word.lower() == groups_to_check[i].lower():
                valid = False
                break
        if not valid:
            continue
        for word in exluded_types: # filters by keyword in name
            if word.lower() in descripts_to_check[i].lower():
                valid = False
                break 
        if valid:
            filtered_cell_types.append(types_to_check[i])

    return filtered_cell_types

def create_df(filepath, columns=None):
    df = pd.read_csv(filepath, sep="\t", header = None)
    if columns == None: # make the first row the column headers
        df.columns = list(df.ix[0])
        df.drop([0], axis=0, inplace=True)
    else:
        df.columns = columns
    return df

def create_dic(df):
    chrms = dict()
    chromosomes = df['Chromosome'].values
    SNPs = df['SNP Name'].values
    for i in range(len(chromosomes)):
        chrm = chromosomes[i]
        if chrm in chrms.keys():
            chrms[chrm].add(SNPs[i])
        else:
            chrms[chrm] = set()
            chrms[chrm].add(SNPs[i])

    allSNPs = []
    allChrms = []
    for chrm in sorted(chrms.keys(), key=get_num):
        print(sorted(chrms[chrm], key=get_num))
        allSNPs.extend(sorted(chrms[chrm], key=get_num))
        allChrms.extend([chrm]*len(chrms[chrm]))
    print(allSNPs)
    counts = dict()
    for i in range(len(SNPs)):
        SNP = SNPs[i]
        if SNP in counts.keys():
            counts[SNP] += 1
        else:
            counts[SNP] = 1

    numCellTypes = []
    for SNP in allSNPs:
        numCellTypes.append(counts[SNP])

    return allSNPs, numCellTypes, allChrms

def plot_enriched_snps(x_vals, y_vals, colors):
    c = {'chr1' :'red',
         'chr2' :'orange',
         'chr3' :'yellow',
         'chr4' :'green',
         'chr5' :'blue',
         'chr6' :'purple',
         'chr7' :'violet',
         'chr8' :'grey',
         'chr9' :'sienna',
         'chr10':'gold',
         'chr11':'olivedrab',
         'chr12':'cyan',
         'chr13':'chocolate',
         'chr14':'fuchsia',
         'chr15':'tan',
         'chr16':'maroon',
         'chr17':'silver',
         'chr18':'salmon',
         'chr19':'dodgerblue',
         'chr20':'hotpink',
         'chr21':'sienna',
         'chr22':'lawngreen'}

    plt.figure(figsize=(20,10))
    plt.barh(list(reversed(range(len(y_vals)))), y_vals, color=[c[x] for x in colors], tick_label=x_vals)
    plt.legend(handles=[mpatches.Patch(color=c[key], label='Chromosome ' + key[3:]) for key in c.keys()], bbox_to_anchor=(1, 1), loc=2)

    plt.ylabel('SNPs in Enriched Regions', fontsize = 15)
    plt.xlabel('Number of Cell Types', fontsize = 15)
    plt.title('SNPs Within Enriched Epigenomic Regions', fontsize = 20)

    #plt.show()
    plt.savefig('./graph_results/SNPs_Within_Enriched_Epigenomic_Regions.png')

def plot_malaria_cell_types(df, meta):
    meta_cell_types = list(meta['Epigenome ID (EID)'].values)

    counts = dict()
    cell_types = df['Cell Type'].values
    filtered_cell_types = filter_cell_types(meta)
    for i in range(len(cell_types)):
        typ = cell_types[i]
        if typ in filtered_cell_types:
            if typ in counts:
                counts[typ] += 1
            else:
                counts[typ] = 1

    x_vals = sorted(counts.keys(), key=get_num) # sort by epigenome name (ie: 'E003')
    #x_vals.sort(key=lambda cell_type: counts[cell_type]) # sort by number of snps
    #x_vals.sort(key=lambda cell_type: meta['GROUP'][meta_cell_types.index(cell_type)]) # sort by group name
    y_vals = [counts[key] for key in x_vals]
    plt.figure(figsize=(15,20))

    plt.barh(list(reversed(range(len(y_vals)))), y_vals, color=['red' if x == 'E066' else meta['COLOR'][meta_cell_types.index(x)] for x in x_vals], tick_label=x_vals)
    
    colors = dict()
    for i in range(3, len(meta_cell_types)):
        if meta['COLOR'][i] in colors.keys():
            continue
        else:
            colors[meta['COLOR'][i]] = meta['GROUP'][i]
    plt.legend(handles=[mpatches.Patch(color='red', label='Liver')] + [mpatches.Patch(color=key, label=colors[key]) for key in colors.keys()], loc=1)

    plt.ylabel('Cell Types', fontsize = 15)
    plt.xlabel('Number of Malaria SNPs', fontsize = 15)
    plt.title('Cell Types With Enriched Malaria SNPs', fontsize = 20)

    #plt.show()
    plt.savefig('./graph_results/Cell_Types_With_Enriched_Malaria_SNPs_[groups].png')

def main():
    df = create_df("./data/epigenomic_overlaps.tsv",['Cell Type', 'Chromosome', 'SNP Start', 'SNP End', 'SNP Name'])
    meta = create_df("./data/roadmap_metadata.tsv")

    #x_vals,y_vals,colors = create_dic(df)
    #plot_enriched_snps(x_vals,y_vals,colors)
    plot_malaria_cell_types(df, meta)

if __name__ == "__main__":
    main()