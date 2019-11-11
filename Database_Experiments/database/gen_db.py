import os
import argparse 
import json
from database import gen_sequences
from database import write_db

def generate(args):
    experiment = args['experiment']
    output_path = args['path']
    output_file_name = args['name']
    window_sizes = args['window_sizes']
    prefix = args['prefix']
    sequences_json = args['sequences_json']

    sequences = None
    with open(sequences_json, 'r') as seqfile:
        sequences = json.load(seqfile)

    seqs = []
    output_files = []
    if 'fractionated' in str(experiment).lower():
        for window_size in window_sizes:
            seqs = []
            seqs = gen_sequences.gen_sequences(sequences['parents']['left_parent']['sequence'], window_size)
            seqs += gen_sequences.gen_sequences(sequences['parents']['right_parent']['sequence'], window_size)
            output_files.append(write_db.write_fasta(output_path, output_file_name + str(window_size), seqs, prefix))

    else: 
        seqs.append(sequences["hybrid"]["sequence"])
        output_files.append(write_db.write_fasta(output_path, output_file_name, seqs, prefix))

    return output_files

