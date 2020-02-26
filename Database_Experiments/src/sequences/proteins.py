import os 
from utils import __file_exists, __make_valid_dir_string, __make_dir
from file_io import fasta
from math import ceil
from random import randint

'''__from_csv

DESC:
    read in proteins from csv into memory
Inputs:
    csv_file: str path to csv file
kwargs:
    has_header: bool if set to true, skip header. Assumed in order name, sequence. Default=False
Outputs:
    list of dictionaries of form {name: str, sequence: str}
'''
def __from_csv(csv_file, has_header=False):
    if not __file_exists(csv_file):
        raise Exception('File {} does not exist'.format(csv_file))
    read_header = False
    prots = []
    with open(csv_file, 'r') as i:
        for line in i:
            if has_header and not read_header:
                read_header = True
                continue
            l = line.split(',')
            prots.append({
                'name': l[0],
                'sequence': l[1]
            })
    return prots

'''load_proteins

DESC:
    return a list of proteins from a source file
Inputs: 
    input_file: str path to a file with proteins. Either csv of form name,sequence or fasta file
Outputs:
    list of dictionaries of the form {'name': str, 'sequence': str}
'''
def load_proteins(input_file):
    if '.csv' in input_file or '.fasta' in input_file:
        prots = __from_csv(input_file) if '.csv' in input_file else fasta.read(input_file) 
    else:
        raise Exception('Protein file should be csv or fasta. File passed in: {}'.format(input_file))
    return prots

'''generate_hybrids

DESC:
    Generate a number of parent proteins with certain constraints
PARMS:
    num_gen: int number of hybrid_proteins to generate
    prots: list of dictionaries of the form
        [{
            'name': str,
            'sequence': str
        }]
kwargs:
    min_contribution: int minimum number of AAs to use from each parent. Default=10
    name_prefix: str prefix to add before every name of protein. Default=HYBRID_
Outputs:
    list of dictionaries. From is
    [{
        left_parent_name: str,
        right_parent_name: str,
        left_parent_sequence: str,
        right_parent_sequence: str,
        left_parent_end: int,
        right_parent_start: int, 
        left_parent_contribution: int,
        right_parent_contribution: int,
        hybrid_protein: str, 
        name: str, 
    }]
'''
def generate_hybrids(prots, num_gen, min_contribution=10, name_prefix='HYBRID_'):
    hybrids = []
    fill_zeros = len(str(ceil(num_gen / 10)))

    for i in range(num_gen):
        print('Generating hybrid protein {}/{}[{}%]\r'.format(i, num_gen, int(float(i)/float(num_gen) * 100)), end="")
        left_parent = prots[randint(0, len(prots)-1)]
        right_parent = prots[randint(0, len(prots)-1)]
        # make sure both parents are long enough
        while len(left_parent['sequence']) < min_contribution:
            left_parent = prots[randint(0, len(prots)-1)]
        while len(right_parent['sequence']) < min_contribution:
            right_parent = prots[randint(0, len(prots)-1)]

        left_end = randint(0, len(left_parent['sequence'])-1)
        right_start = randint(0, len(right_parent['sequence'])-1)
        while left_end < min_contribution:
            left_end = randint(0, len(left_parent['sequence'])-1)
        while len(right_parent['sequence']) - right_start < min_contribution:
            right_start = randint(0, len(right_parent['sequence'])-1)
        hybrid = left_parent['sequence'][:left_end] + right_parent['sequence'][right_start:]
        name = name_prefix + str(i).zfill(fill_zeros)
        hybrid_d = {
            'left_parent_name': left_parent['name'],
            'right_parent_name': right_parent['name'],
            'left_parent_sequence': left_parent['sequence'],
            'right_parent_sequence': right_parent['sequence'],
            'left_parent_end': left_end,
            'right_parent_start': right_start,
            'left_parent_contribution': left_end - 1,
            'right_parent_contribution': len(right_parent['sequence']) - right_start,
            'protein': hybrid,
            'name': name
        }
        hybrids.append(hybrid_d)
    print('\nFinished generating hybrid proteins')
    return hybrids

'''generate

DESC:
    generates fasta files for k-mers
Inputs:
    peptides: list of dictionaries of form {'name': str, 'sequence': str}
kwargs:
    save_dir: str file path to directory to save to. Default=./
Outputs:
    list of strings of files created
'''
def generate_databases(peptides, save_dir='./'):
    output_files = []
    save_dir = __make_valid_dir_string(save_dir) + 'databases/'
    __make_dir(save_dir)

    for pep in peptides:
        file_name = '{}{}'.format(save_dir, pep['name'])
        output_files.append(fasta.write(file_name, [pep]))
    return output_files
