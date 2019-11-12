import argparse 
import json
from spectra import gen_spectra, write_spectra
import sys 
sys.path.append('../')
from database import gen_sequences


def generate(args):
    experiment = args['experiment']
    output_path = args['path']
    output_file_name = args['name']
    window_sizes = args['window_sizes']
    title_prefix = args['title_prefix']
    sequences_json = args['sequences_json']

    sequences = None
    with open(sequences_json, 'r') as seqfile:
        sequences = json.load(seqfile)

    seqs = []
    output_files = []

    if 'fractionated' in str(experiment).lower():
        seqs.append(sequences["hybrid"]["sequence"])
        spectra = gen_spectra.gen_spectra(seqs)
        output_files.append(write_spectra.write_mzml(output_path, output_file_name, spectra, title_prefix))

    else: 
        for window_size in window_sizes:
            seqs = gen_sequences.gen_sequences(sequences['parents']['left_parent']['sequence'], window_size)
            seqs += gen_sequences.gen_sequences(sequences['parents']['right_parent']['sequence'], window_size)
            spectra = gen_spectra.gen_spectra(seqs)
            output_files.append(write_spectra.write_mzml(output_path, output_file_name + str(window_size), spectra, title_prefix))
    return output_files

# if __name__ == '__main__':
#     parser = argparse.ArgumentParser(description='generate a database and store in a file with which experiment')
#     parser.add_argument('experiment', metavar='E', type=str, help='Which experiment to generate a database for. Options are Flipped or Fractionated Defaults to Flipped.')
#     parser.add_argument('path', metavar='P', type=str, help='Path in which to save the database')
#     parser.add_argument('name', metavar='N', type=str, help='Name to give to the new database')
#     # parser.add_argument('output_format', metavar='O', type=str, help='Output format. Only current one is mgf')
#     parser.add_argument('--title_prefix', dest='title_prefix', default='Spectrum ', help='Prefix for the title of the spectrum')
#     parser.add_argument('--window_size', dest='window_size', default=3, help='Size of the window to use in generating the fractionated db. Defaults to 3')
#     args = parser.parse_args()

#     gen(args)
