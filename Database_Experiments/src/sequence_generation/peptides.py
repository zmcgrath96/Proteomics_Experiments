from sequence_generation import digest
from utils import __file_exists

'''gen_peptides

DESC:
    generate peptides from proteins
PARAMS:
    prots: list of dictionaries of the form {'name': str, 'sequence': str}
OPTIONAL:
    number_peptides: int number of peptides to generate. Default=10
    min_length: int minimum length peptide to generate. Default=3
    max_length: int maximum length peptide to generate. Default=20
    save_dir: str path to file to save peptides. Default=./ 
RETURNS:
    list of dictionaries of form 
    {
        'peptide_name': str,
        'peptide_sequence': str,
        'parent_name': str,
        'parent_sequence': str,
        'start_index': int, 
        'end_index': int
    }
'''
def gen_peptides(prots, number_peptides=10, min_length=3, max_length=20, save_dir='./'):
    return digest.tryptic(prots, number_peptides, min_length=min_length, max_length=max_length, save_dir=save_dir)
    
'''load_peptides

DESC:
    load peptides from file
PARAMS:
    peptide_file: str file path to peptide file
RETURNS:
    list of dictionaries of form 
    {
        'peptide_name': str,
        'peptide_sequence': str,
        'parent_name': str,
        'parent_sequence': str,
        'start_index': int, 
        'end_index': int
    }
'''
def load_peptides(peptide_file):
    if not __file_exists(peptide_file):
        Exception('Digestion file should be a valid file. {} is not a file'.format(peptide_file))
    return digest.load_digest(peptide_file)

'''get_peptides

DESC: 
    bring peptides into memory
OPTIONAL:
    prots: list of dictionaries of the form {'name': str, 'sequence': str}. Default=[]
    number_peptides: int number of peptides to generate. Default=10
    min_length: int minimum length peptide to generate. Default=3
    max_length: int maximum length peptide to generate. Default=20
    save_dir: str path to file to save peptides. Default=./ 
    peptide_file:  str file path to peptide file. Default=''
RETURNS: 
    list of dictionaries of form 
    {
        'peptide_name': str,
        'peptide_sequence': str,
        'parent_name': str,
        'parent_sequence': str,
        'start_index': int, 
        'end_index': int
    }
'''
def get_peptides(prots=[], number_peptides=10, min_length=3, max_length=20, save_dir='./', peptide_file=''):
    if not peptide_file or peptide_file == '' or not __file_exists(peptide_file):
        return gen_peptides(prots, number_peptides=number_peptides, min_length=min_length, max_length=max_length, save_dir=save_dir)
    else:
        return load_peptides(peptide_file) 
