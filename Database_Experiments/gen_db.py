import os
import argparse 
import json
import gen_sequences
import write_fasta

def gen(args):
    experiment = args.experiment
    output_path = args.path
    output_file_name = args.name
    window_size = int(args.window_size)
    prefix = args.prefix

    sequences = None
    with open('./sequences.json', 'r') as seqfile:
        sequences = json.load(seqfile)

    seqs = []
    if experiment == 'Fractionated':
        seqs = gen_sequences.gen_sequences(sequences['parents']['left_parent']['sequence'], window_size)
        seqs += gen_sequences.gen_sequences(sequences['parents']['right_parent']['sequence'], window_size)

    else: 
        seqs.append(sequences["hybrid"]["sequence"])

    write_fasta.write(output_path, output_file_name, seqs, prefix)
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='generate a database and store in a file with which experiment')
    parser.add_argument('experiment', metavar='E', type=str, help='Which experiment to generate a database for. Options are Flipped or Fractionated Defaults to Flipped.')
    parser.add_argument('path', metavar='P', type=str, help='Path in which to save the database')
    parser.add_argument('name', metavar='N', type=str, help='Name to give to the new database')
    parser.add_argument('--window_size', dest='window_size', default=3, help='Size of the window to use in generating the fractionated db. Defaults to 3')
    parser.add_argument('--prefix', dest='prefix', default='HYBRID_', help='Prefix to add to name of the sequences')
    args = parser.parse_args()

    gen(args)
