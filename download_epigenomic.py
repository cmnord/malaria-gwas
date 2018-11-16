"""
Docstring
"""
import requests
import os.path
import pdb

BASE_URL = 'https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/{}'
FILE_NAME = 'epigenomic_annotations/{}-H3K27ac.narrowPeak.gz'


def get_cell_types():
    """
    Docstring
    """
    cell_types = []
    with open('epigenomic_annotations/cell_types.txt', 'r') as cell_file:
        for line in cell_file:
            # remove newline
            cell_types.append(line[:-1])

    return cell_types


def main():
    """
    Docstring
    """
    cell_types = ['E{:03}'.format(i) for i in range(3, 130)]
    for cell_type in cell_types:
        filename = FILE_NAME.format(cell_type)
        print(cell_type)
        if os.path.isfile(filename):
            print('skip')
            continue
        url = BASE_URL.format(filename)
        r = requests.get(url)
        if r.status_code == '200':
            open(filename, 'wb').write(r.content)


if __name__ == "__main__":
    main()
