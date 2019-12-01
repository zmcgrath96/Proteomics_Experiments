import os
import argparse 
import json
from src.database import gen_sequences, write_db

def generate(args):
    experiment = args['experiment']
    output_path = args['path']
    output_file_name = args['name']
    window_sizes = args['window_sizes']
    prefix = args['prefix']
    sequences = args['sequences_dict']
    peptide_index = args['peptide_index']

    seqs = []
    output_files = []
    if 'fractionated' in str(experiment).lower():
        for window_size in window_sizes:
            seqs = []
            parent_one_output_name = '{}_{}_{}'.format(output_file_name, sequences['parents']['left_parent']['name'], window_size)
            parent_two_output_name = '{}_{}_{}'.format(output_file_name, sequences['parents']['right_parent']['name'], window_size)
            seqs = gen_sequences.gen_sequences(sequences['parents']['left_parent']['sequence'], window_size)
            output_files.append(write_db.write_fasta(output_path, parent_one_output_name, seqs, prefix))
            seqs = gen_sequences.gen_sequences(sequences['parents']['right_parent']['sequence'], window_size)
            output_files.append(write_db.write_fasta(output_path, parent_two_output_name, seqs, prefix))

    else: 
        seqs += [sequences['hybrid']['sequence']]
        output_files.append(write_db.write_fasta(output_path, output_file_name, seqs, prefix))
        seqs = sequences[peptide_index]
        output_name = 'peptide_{}'
        for i, seq in enumerate(seqs):
            output_files.append(write_db.write_fasta(output_path, output_name.format(i), [seq], prefix))

    return output_files

