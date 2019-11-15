# Generate target sequences
from pyteomics import fasta
import random
import os
import json

digests = None

def load_fasta(db_file_path):
    entries = []
    with fasta.read(db_file_path) as db:
        for entry in db:
            entries.append(entry)
    return entries

def load_digests():
    cwd = os.path.dirname(os.path.realpath(__file__))
    parent_dir_idx = cwd.split('/').index('Proteomics_Experiments')
    digest_json = '/'.join(cwd.split('/')[:parent_dir_idx].append('utils/digests.json'))
    digest = None
    with open(digest_json, 'r') as o:
        digest = json.load(o)
    return digest

def digest(protein_seq, digest, missed_cleavages, missed_cleavage_prob):
    if digests is None:
        digests = load_digests()
    random_pos =  random.randint(0, len(protein_seq) - 1)
    digest_dict = {}
    for x in digests[digest]:
        digest_dict[x['site']] = x['terminus']
    # see if were on a cleaving site
    if protein_seq[random_pos] in digest_dict:
        #pick left or right 
        if random.random() > .5: 
            # go right
            pass 
        else:
            # go left
            pass
    else:
        #go left and right until we hit a cleaving site
        pass

'''generate

DESC:
    Given a reference database generate some possible sequences 
    Uses the digest to determine cleaving sites 
PARAMS:
    ref_db: string filepath to the reference database
    number: int number of sequences to generate
OPTIONAL:
    digest: string the digest to perfrom on the proteins. Defaults to trypsin
    missed_cleavages: int number of cleavages to potentially miss
    missed_cleavage_prob: float number between 0, 1 to determine how often to miss a cleaving site
RETRUNS: 
    list of tuples with the generated peptide and the parent protein (peptide, parent protein sequence, parent protein name)
'''
def generate(ref_db, number, digest='trypsin', missed_cleavages=0, missed_cleavage_prob=0):
    proteins = load_fasta(ref_db)
    peptides = []
    for _ in range(number):
        random_protein_indx = random.randint(0, len(proteins))
        parent_sequence = proteins[random_protein_indx].sequence
        parent_desc = proteins[random_protein_indx].description
        peptide = digest(parent_sequence, digest, missed_cleavages, missed_cleavage_prob)
        tup = (peptide, parent_sequence, parent_desc)
        peptides.append(tup)
    return peptides

import sys 
generate(sys.argv[1], int(sys.argv[2]))