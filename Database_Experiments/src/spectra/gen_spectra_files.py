import argparse 
import json
from spectra import gen_spectra, write_spectra
from sequence_generation.gen_k_mers import k_mers


def generate(args):
    output_path = args['path']
    output_file_name = args['name']
    window_sizes = args['window_sizes']
    title_prefix = args['title_prefix']
    sequences = args['sequences_dict']

    seqs = []
    output_files = []
    
    for window_size in window_sizes:
        print('Generating {}-mer spectra for all proteins...'.format(window_size))
        for sequence in sequences['sample']['proteins']:
            name = '{}_{}_{}'.format(output_file_name, sequence['name'], window_size)
            seqs = gen_sequences.gen_sequences(sequence['sequence'], window_size)
            spectra = gen_spectra.gen_spectra(seqs)
            output_files.append(write_spectra.write_mzml(output_path, name, spectra, title_prefix))

    return output_files
