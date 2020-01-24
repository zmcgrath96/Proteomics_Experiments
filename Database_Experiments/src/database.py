import os
import argparse 
import json
from utils import __make_valid_dir_string, __make_dir
from file_io import fasta
from math import ceil

'''generate

DESC:
    generates fasta files for k-mers
PARAMS:
    peptides: list of dictionaries of form {'name': str, 'sequence': str}
OPTIONAL:
    save_dir: str file path to directory to save to. Default=./
RETURNS:
    list of strings of files created
'''
def generate(peptides, save_dir='./'):
    output_files = []
    save_dir = __make_valid_dir_string(save_dir) + 'databases/'
    __make_dir(save_dir)

    for pep in peptides:
        file_name = '{}{}'.format(save_dir, pep['name'])
        output_files.append(fasta.write(file_name, [pep]))
    return output_files

