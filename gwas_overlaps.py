"""
Find overlaps between:
1. Significant SNPs from malaria GWAS and
    A. The epigenomic state of different cell types (specifically H3K4me3)
    B. The gene expression level of red blood cells
2. Background SNPs and
    A. The epigenomic state of different cell types (specifically H3K4me3)
    B. The gene expression level of red blood cells
"""
import subprocess
import os

MALARIA_GWAS = 'data/gwas_associations/snps.bed'
BACKGROUND_GWAS = 'data/gwas_associations/gwas_catalog_data.bed'
RBC_EXPRESSION_FILE = 'data/rbc_gene_expression/rbc_w_genes_locs.bed'
EPIGENOMIC_FOLDER = 'data/epigenomic_annotations/H3K4me3/'
OUTPUT_FOLDER = 'data/overlaps/H3K4me3/{}'
ZIPFILE_EXTENSION = '.gz'


def is_zipped(filename):
    return filename[-len(ZIPFILE_EXTENSION):] == ZIPFILE_EXTENSION


def get_cell_type(filename):
    # filename format is E003-H3K27ac.narrowPeak
    return filename[:4]


def acetylation_overlap(gwas_file, output_file):
    epigenomic_files = sorted(os.listdir(EPIGENOMIC_FOLDER))
    command = ['bedtools', 'intersect', '-a', gwas_file, '-b']
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


def rbc_expression_overlap(gwas_file, output_file):
    with open(output_file, 'w') as f:
        out = subprocess.run(
            ['bedtools', 'intersect', '-a', RBC_EXPRESSION_FILE, '-b', gwas_file], stdout=subprocess.PIPE
        )
        results = out.stdout.decode('utf-8').strip()
        f.write(results)


if __name__ == "__main__":
    acetylation_overlap(MALARIA_GWAS, OUTPUT_FOLDER.format(
        'malaria_acetylation_overlap.tsv'))
    rbc_expression_overlap(
        MALARIA_GWAS, OUTPUT_FOLDER.format('malaria_rbc_expression_overlap.tsv'))
    acetylation_overlap(
        BACKGROUND_GWAS, OUTPUT_FOLDER.format('background_acetylation_overlap.tsv'))
    rbc_expression_overlap(
        BACKGROUND_GWAS, OUTPUT_FOLDER.format('background_rbc_expression_overlap.tsv'))
