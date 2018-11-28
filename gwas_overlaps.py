"""
Find overlaps between significant SNPs from malaria GWAS and
1. The epigenomic state of different cell types (specifically H3K27 acetylation)
2. The gene expression level of red blood cells
"""
import subprocess
import os
import pdb

GWAS_FILE = 'data/gwas_associations/snps.bed'
ZIPFILE_EXTENSION = '.gz'


def is_zipped(filename):
    return filename[-len(ZIPFILE_EXTENSION):] == ZIPFILE_EXTENSION


def get_cell_type(filename):
    # filename format is E003-H3K27ac.narrowPeak
    return filename[:4]


def calculate_acetylation_overlap():
    EPIGENOMIC_FOLDER = 'data/epigenomic_annotations/'
    epigenomic_files = sorted(os.listdir(EPIGENOMIC_FOLDER))
    command = ['bedtools', 'intersect', '-a', GWAS_FILE, '-b']
    output_file = 'data/epigenomic_overlaps.tsv'
    with open(output_file, 'w') as f:
        for filename in epigenomic_files:
            if is_zipped(filename):
                continue

            cell_type = get_cell_type(filename)
            print(cell_type)

            out = subprocess.run(
                command + [EPIGENOMIC_FOLDER+filename], stdout=subprocess.PIPE)

            if len(out.stdout) == 0:
                continue

            results = out.stdout.decode('utf-8').strip().split('\n')
            for line in results:
                f.write('{}\t{}\n'.format(cell_type, line))


def calculate_rbc_expression_overlap():
    output_file = 'data/rbc_expression_overlaps.tsv'
    RBC_EXPRESSION_FILE = 'data/rbc_gene_expression/rbc_w_genes_locs.bed'
    with open(output_file, 'w') as f:
        out = subprocess.run(
            ['bedtools', 'intersect', '-a', RBC_EXPRESSION_FILE, '-b', GWAS_FILE], stdout=subprocess.PIPE
        )
        results = out.stdout.decode('utf-8').strip()
        f.write(results)


if __name__ == "__main__":
    # calculate_acetylation_overlap()
    calculate_rbc_expression_overlap()
