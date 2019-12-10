from sequence_generation import digest
from utils import __file_exists

def gen_peptides(sequence_dict, number_peptides, peptide_index='peptides', min_length=3, max_length=20, save_dir='./'):
    proteins = sequence_dict['sample']['proteins']
    peptides = digest.tryptic(proteins, number_peptides, min_length=min_length, max_length=max_length, save_dir=save_dir)
    sequence_dict[peptide_index] = peptides

def load_peptides(sequence_dict, digestion_file, peptide_index='peptides'):
    if not __file_exists(digestion_file):
        Exception('Digestion file should be a valid file. {} is not a file'.format(digestion_file))
    d = digest.load_digest(digestion_file)
    peps = [x['peptide_sequence'] for _, x in d.items()]
    sequence_dict[peptide_index] = peps
