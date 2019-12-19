import os
import argparse 
import json
from sequence_generation.gen_k_mers import k_mers
from utils import __make_valid_dir_string, __make_dir, __make_valid_fasta_file
from math import ceil

############################################################
#                   "PRIVATE" FUNCTIONS
############################################################
'''__write_fasta

DESC:
    write a fasta file
PARAMS:
    output_name: str name of file to write to 
    sequences: list of dictionaries of form {'name': str, 'sequence': str}
RETURNS:
    name of the output file written to
'''
def __write_fasta(output_name, sequences):
    output_name = __make_valid_fasta_file(output_name)
    with open(output_name, 'w') as o:
        for seq in sequences:
            o.write('>{}\n{}\n'.format(seq['name'], seq['sequence']))
    return output_name

############################################################
#                END "PRIVATE" FUNCTIONS
############################################################

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
        output_files.append(__write_fasta(file_name, [pep]))
    return output_files

