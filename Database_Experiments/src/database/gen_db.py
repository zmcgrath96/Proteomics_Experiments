import os
import argparse 
import json
from database import gen_sequences, write_db

def generate(args):
    output_path = args['path']
    output_file_name = args['name']
    prefix = args['prefix']
    sequences = args['sequences_dict']
    peptide_index = args['peptide_index']

    seqs = []
    output_files = []

    seqs += [sequences['hybrid']['sequence']]
    output_files.append(write_db.write_fasta(output_path, output_file_name, seqs, prefix))
    seqs = sequences[peptide_index]
    output_name = 'peptide_{}'
    for i, seq in enumerate(seqs):
        num = str(i) if i > 9 else '0' + str(i)
        output_files.append(write_db.write_fasta(output_path, output_name.format(num), [seq], prefix))

    return output_files

