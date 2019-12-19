import argparse 
import json
from spectra import gen_spectra, write_spectra
from sequence_generation.gen_k_mers import k_mers
from utils import __make_valid_dir_string, __make_dir

'''generate 

DESC:
    generate kmers and the spectra files associated with them 
PARAMS:
    sequences: list of dictionaries of the form {'name': str, 'sequence': str}. Sequences to generate
                kmers from
    window_sizes: list of ints size of kmers to generate
OPTIONAL:
    save_dir: str the directory in which to save all the spectra files. Default=./
RETURNS:
    list of strs of the file names/paths generated
'''
def generate(sequences, window_sizes, save_dir='./'):
    output_files = []
    save_dir = __make_valid_dir_string(save_dir) + 'spectra/'
    
    for window_size in window_sizes:
        print('Generating {}-mer spectra for all proteins...'.format(window_size))
        for sequence in sequences:
            name = '{}_{}'.format(sequence['name'], window_size)
            kmers = k_mers(sequence['sequence'], window_size)
            spectra = gen_spectra.gen_spectra(kmers)
            output_files.append(write_spectra.write_mzml(name, spectra, output_dir=save_dir))

    return output_files
