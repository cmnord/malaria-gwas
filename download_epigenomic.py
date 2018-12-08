"""
Docstring
"""
import requests
import os.path
import pdb

BASE_URL = 'https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/{}'
FILE_NAME = '{}-H3K4me3.narrowPeak.gz'


def main():
    """
    Docstring
    """
    cell_types = ['E{:03}'.format(i) for i in range(1, 130)]
    for cell_type in cell_types:
        filename = FILE_NAME.format(cell_type)
        print(filename)
        print(cell_type)
        if os.path.isfile(filename):
            print('skip')
            continue
        url = BASE_URL.format(filename)
        r = requests.get(url)
        if r.status_code != 200:
            print(r.status_code)
            continue
        open(filename, 'wb').write(r.content)


if __name__ == "__main__":
    main()
