import sys
import pandas as pd
import matplotlib.patches as mpatches
from matplotlib import pyplot as plt


def get_num(elem):
    for i in range(len(elem)):
        if elem[i] == 'X':
            return -1
        elif elem[i] == 'Y':
            return 0
        try:
            int(elem[i])
            return int(elem[i:])
        except:
            continue


def add_SNP_names(df, SNPs):
    df.rename(columns={'SNP Name': 'GeneID'}, inplace=True)

    names = []
    for i in range(len(df['SNP Start'])):
        names.append(
            SNPs['SNP'][[int(x) for x in SNPs['Position'].values].index(df['SNP Start'][i])])

    df['SNP Name'] = names


def filter_cell_types(meta):
    filtered_meta = meta.copy()

    exluded_groups = {'iPSC'}
    exluded_types = {'fetal', 'culture', 'carcinoma', 'leukemia'}

    filtered_cell_types = []
    types_to_check = meta['Epigenome ID (EID)'].values[2:]
    groups_to_check = meta['GROUP'].values[2:]
    descripts_to_check = meta['Standardized Epigenome name'].values[2:]

    for i in range(len(types_to_check)):
        valid = True
        for word in exluded_groups:  # filters by group
            if word.lower() == groups_to_check[i].lower():
                valid = False
                filtered_meta.drop([i+3], axis=0, inplace=True)
                break
        if not valid:
            continue
        for word in exluded_types:  # filters by keyword in name
            if word.lower() in descripts_to_check[i].lower():
                valid = False
                filtered_meta.drop([i+3], axis=0, inplace=True)
                break
        if valid:
            filtered_cell_types.append(types_to_check[i])

    # creates new metadata that's filtered for valid cell types
    filtered_meta.to_csv(
        "./data/roadmap_metadata_filtered.tsv", "\t", index=False)

    return filtered_cell_types


def create_df(filepath, columns=None):
    if filepath[-3:] == "tsv":
        df = pd.read_csv(filepath, sep="\t", header=None)
    elif filepath[-3:] == "csv":
        df = pd.read_csv(filepath, sep=",", header=None)
    if columns == None:  # make the first row the column headers
        df.columns = list(df.ix[0])
        df.drop([0], axis=0, inplace=True)
    else:
        df.columns = columns
    return df


def create_dic(primary, blood):
    chrms1 = dict()
    chromosomes = primary['Chromosome'].values
    primeSNPs = primary['SNP Name'].values
    for i in range(len(chromosomes)):
        chrm = chromosomes[i]
        if chrm in chrms1.keys():
            chrms1[chrm].add(primeSNPs[i])
        else:
            chrms1[chrm] = set()
            chrms1[chrm].add(primeSNPs[i])

    chrms2 = dict()
    chromosomes = blood['Chromosome'].values
    bloodSNPs = blood['SNP Name'].values
    for i in range(len(chromosomes)):
        chrm = chromosomes[i]
        if chrm in chrms2.keys():
            chrms2[chrm].add(bloodSNPs[i])
        else:
            chrms2[chrm] = set()
            chrms2[chrm].add(bloodSNPs[i])

    allSNPs = []
    allChrms = []
    for chrm in sorted(set(list(chrms1.keys()) + list(chrms2.keys())), key=get_num):
        if chrm in chrms1.keys():
            allSNPs.extend(sorted(chrms1[chrm], key=get_num))
            allChrms.extend([chrm]*len(chrms1[chrm]))
        if chrm in chrms2.keys():
            allSNPs.extend(sorted(chrms2[chrm]))
            allChrms.extend([chrm]*len(chrms2[chrm]))

    counts = dict()
    for i in range(len(primeSNPs)):
        SNP = primeSNPs[i]
        if SNP in counts.keys():
            counts[SNP] += 1
        else:
            counts[SNP] = 1
    for i in range(len(bloodSNPs)):
        SNP = bloodSNPs[i]
        if SNP in counts.keys():
            counts[SNP] += 1
        else:
            counts[SNP] = 1

    numCellTypes = []
    for SNP in allSNPs:
        numCellTypes.append(counts[SNP])

    return allSNPs, numCellTypes, allChrms


def plot_enriched_snps(x_vals, y_vals, colors):
    c = {'chr1': 'red',
         'chr2': 'orange',
         'chr3': 'yellow',
         'chr4': 'green',
         'chr5': 'blue',
         'chr6': 'purple',
         'chr7': 'violet',
         'chr8': 'grey',
         'chr9': 'sienna',
         'chr10': 'gold',
         'chr11': 'olivedrab',
         'chr12': 'cyan',
         'chr13': 'chocolate',
         'chr14': 'fuchsia',
         'chr15': 'tan',
         'chr16': 'maroon',
         'chr17': 'silver',
         'chr18': 'salmon',
         'chr19': 'dodgerblue',
         'chr20': 'hotpink',
         'chr21': 'sienna',
         'chr22': 'lawngreen',
         'chrX': 'lightcoral',
         'chrY': 'crimson'}

    plt.figure(figsize=(20, 10))
    plt.barh(list(reversed(range(len(y_vals)))), y_vals,
             color=[c[x] for x in colors], tick_label=x_vals)
    plt.legend(handles=[mpatches.Patch(color=c[key], label='Chromosome ' + key[3:])
                        for key in c.keys()], bbox_to_anchor=(1, 1), loc=2)

    for i, v in enumerate(reversed(y_vals)):
        plt.text(v, i, " " + str(v), color='black',
                 va='center', fontweight='bold')

    plt.ylabel('SNPs in Enriched Regions', fontsize=15)
    plt.xlabel('Number of Cell Types', fontsize=15)
    plt.title('SNPs Within Enriched Epigenomic Regions', fontsize=20)

    # plt.show()
    plt.savefig(
        './graph_results/All_SNPs_Within_Enriched_Epigenomic_Regions.png')


def plot_malaria_cell_types(primary, background, blood, blood_back, meta):
    meta_cell_types = ['E000'] + list(meta['Epigenome ID (EID)'].values)

    counts1 = dict()
    cell_types = primary['Cell Type'].values
    filtered_cell_types = filter_cell_types(meta)
    for i in range(len(cell_types)):
        typ = cell_types[i]
        if typ in filtered_cell_types:
            if typ in counts1:
                counts1[typ] += 1
            else:
                counts1[typ] = 1
    counts1['E000'] = len(blood['SNP Name'].values)

    counts2 = dict()
    cell_types = background['Cell Type'].values
    for i in range(len(cell_types)):
        typ = cell_types[i]
        if typ in filtered_cell_types:
            if typ in counts2:
                counts2[typ] += 1
            else:
                counts2[typ] = 1
    counts2['E000'] = len(blood_back['SNP Name'].values)

    mix = dict()
    for key in counts1.keys():
        mix[key] = float(counts1[key])/float(counts1[key]+counts2[key])

    # sort by epigenome name (ie: 'E003')
    x_vals = sorted(counts1.keys(), key=get_num)
    # x_vals.sort(key=lambda cell_type: mix[cell_type]) # sort by number of snps
    x_vals.sort(key=lambda cell_type: (['Erythrocytes (RBC)'] + list(
        meta['GROUP'].values))[meta_cell_types.index(cell_type)])  # sort by group name
    y_vals = [mix[key] for key in x_vals]
    plt.figure(figsize=(15, 20))

    plt.barh(list(reversed(range(len(y_vals)))), y_vals, color=[
             'blue' if x == 'E066' else 'red' if x == "E000" else meta['COLOR'][meta_cell_types.index(x)] for x in x_vals], tick_label=x_vals)

    colors = dict()
    for i in range(3, len(meta_cell_types)):
        if meta['COLOR'][i] in colors.keys():
            continue
        else:
            colors[meta['COLOR'][i]] = meta['GROUP'][i]
    plt.legend(handles=[mpatches.Patch(color='red', label='Erythrocytes (RBC)'), mpatches.Patch(
        color='blue', label='Liver')] + [mpatches.Patch(color=key, label=colors[key]) for key in colors.keys()], loc=1)

    # for i, v in enumerate(reversed(y_vals)):
    #     plt.text(v, i, " " + str(v), color='black', va='center', fontweight='bold')

    plt.ylabel('Cell Types', fontsize=15)
    plt.xlabel('Percent Malaria SNPs of All SNPs in Enriched Regions', fontsize=15)
    plt.title('Cell Types With Enriched Malaria SNPs', fontsize=20)

    # plt.show()
    plt.savefig(
        './graph_results/All_Cell_Types_With_Enriched_Malaria_SNPs_[normalized].png')


def main():
    primary = create_df("./data/overlaps/malaria_acetylation_overlap.tsv",
                        ['Cell Type', 'Chromosome', 'SNP Start', 'SNP End', 'SNP Name'])
    back = create_df("./data/overlaps/background_acetylation_overlap.tsv",
                     ['Cell Type', 'Chromosome', 'SNP Start', 'SNP End', 'SNP Name'])
    p_meta = create_df("./data/roadmap_metadata.tsv")

    blood = create_df("./data/overlaps/malaria_rbc_expression_overlap.tsv",
                      ['Chromosome', 'SNP Start', 'SNP End', 'SNP Name'])
    add_SNP_names(blood, create_df("./data/snps.csv"))
    blood_back = create_df("./data/overlaps/background_rbc_expression_overlap.tsv",
                           ['Chromosome', 'SNP Start', 'SNP End', 'SNP Name'])

    # finding how many cell types each malaria SNP is enriched in
    x_vals, y_vals, colors = create_dic(primary, blood)
    plot_enriched_snps(x_vals, y_vals, colors)

    # finding how many enriched malaria SNPs (out of all enriched SNPs) are in each cell type
    plot_malaria_cell_types(primary, back, blood, blood_back, p_meta)


if __name__ == "__main__":
    main()
