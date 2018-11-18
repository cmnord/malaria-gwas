"""
Find overlaps between significant SNPs from malaria GWAS
and the epigenomic state of different cell types
(specifically H3K27 acetylation).
"""
import subprocess
import os
import pdb

gwas_file = 'data/gwas_associations/snps.bed'
epigenomic_folder = 'data/epigenomic_annotations/'
output_file = 'data/epigenomic_overlaps.tsv'

zipped_end = '.gz'
def is_zipped(filename):
    return filename[-len(zipped_end):] == zipped_end

def get_cell_type(filename):
    # filename format is E003-H3K27ac.narrowPeak
    return filename[:4]

def make_table():
    command = ['bedtools', 'intersect', '-a', gwas_file, '-b']
    epigenomic_files = sorted(os.listdir(epigenomic_folder))
    with open(output_file, 'w') as f:
        for filename in epigenomic_files:
            if is_zipped(filename):
                continue

            cell_type = get_cell_type(filename)
            print(cell_type)

            out = subprocess.run(command + [epigenomic_folder+filename], stdout=subprocess.PIPE)

            if len(out.stdout) == 0:
                continue

            results = out.stdout.decode('utf-8').strip().split('\n')
            for line in results:
                f.write('{}\t{}\n'.format(cell_type, line))


if __name__ == "__main__":
    make_table()