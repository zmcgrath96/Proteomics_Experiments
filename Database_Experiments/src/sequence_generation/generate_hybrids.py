from utils import __file_exists
from sequence_generation import hybrids

'''generatate hybrids

DESC:
    bring hybrid proteins and hybrid peptides into memory
PARAMS:
OPTIONAL:
    hyb_pep_file: str path to a hybrid peptide file. If none is given, new ones generated. Default=''
    hyb_prot_file: str path to a hybrid protein file. If none is given, new ones generated. Default=''
    prots: list of dictionaries of form {'name': str, 'sequence': str}. Used if new ones have to be generated. Default=[]
    num_gen: int number of hybrids to generate. Default=10
    min_contribution: int minimum number from each protein to take for generation. Default=10
    min_length: int minimum length peptide to generate. Default=10
    max_length: int maximum length peptide to generate. Default=20
    save_dir: str path to directory in which to save the hybrid proteins and peptides. Default=./
RETURNS:
    l1, l2:
        l1: list of hybrid peptides
        l2: list of hybrid proteins
'''
def generate_hybrids(hyb_pep_file='', hyb_prot_file='', prots=[], num_gen=10, min_contribution=10, min_length=10, max_length=20, save_dir='./'):
    hyb_prots = None 
    hyb_peps = None 
    if not __file_exists(hyb_prot_file):
        hyb_prots = hybrids.generate_proteins(num_gen, prots, min_contribution=10)
    else:
        hyb_prots = hybrids.read_proteins(hyb_prot_file)
    
    if not __file_exists(hyb_pep_file):
        hyb_peps, _ = hybrids.generate_peptides(hyb_prots, num_gen=num_gen, min_length=min_length, max_length=max_length)
    else:
        hyb_peps, _ = hybrids.read_peptides(hyb_pep_file)

    return hyb_peps, hyb_prots
