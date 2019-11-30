from random import randint, choice, random
import sys 
sys.path.append('../')
from utils import __make_dir, __make_valid_dir_string

'''load_digest

DESC:
    load the digests from the digestion tsv
PARAMS:
    digest_file: string file path for the digestion tsv
RETURNS:
    dictionary: keys are peptide names and entries have 
    {
        'peptide_sequence': string,
        'parent_name': string,
        'parent_sequence': string,
        'start_index': int
    }
'''
def load_digest(digest_file):
    digests = {}
    with open(digest_file, 'r') as o:
        for i, line in enumerate(o):
            if i == 0:
                continue
            l = line.split('\t')
            name = 'peptide_' + str(i-1)
            entry = {'peptide_sequence': l[0], 'parent_name': l[1], 'parent_sequence': l[2], 'start_index': int(l[3])}
            digests[name] = entry 
    return digests

'''__tryptic_digest

DESC:
    Perfrom the actual digestion 
PARAMS: 
    sequence: full length sequence to perform digestion on
'''
def __tryptic_digest(sequence, miss_prob):
    # pick a start point thats not the start or the end
    start = randint(0, len(sequence)-1)
    while start == 0 or start == len(sequence) - 1:
        start = randint(0, len(sequence) - 1)
    end = start

    # if start is at a sequence point, go left or right
    if sequence[start] == 'R' or sequence[start] == 'K':
        # go left
        if random() < 0.5:
            while start > 0:
                start -= 1
                if sequence[start] == 'K' or sequence[start] == 'R':
                    if random() < miss_prob:
                        continue
                    break 
            this_pep = sequence[start+1:end+1]
        # go right
        else:
            while end < len(sequence) - 1:
                end += 1
                if sequence[end] == 'K' or sequence[end] == 'R':
                    if random() < miss_prob:
                        continue
                    break 
            this_pep = sequence[start+1:end+1]

    # go both ways
    else:
        while start > 0:
            if sequence[start] == 'K' or sequence[start] == 'R':
                if random() < miss_prob:
                    continue
                break 
            start -= 1
        while end < len(sequence) - 1:
            if sequence[end] == 'K' or sequence[end] == 'R':
                if random() < miss_prob:
                    continue
                break 
            end += 1
        this_pep = sequence[start+1:end+1]
    
    return this_pep, start + 1

'''tryptic

DESC:
    Generate peptides with tryptic digest from sequences with missed cleavages at a probability
PARAMS:
    sequences: list of dictionary objects with entries {'sequences': string, 'name': string}
    number_digests: int number of peptides to generate
OPTIONAL:
    miss_prob: float the probability of a missed cleavage. Should range [0, 1). Default=0
    save_dir: string the directory in which to save the digestion file. Default=./
    save_name: string the name to save the digestion information in. Default=digestion.tsv
    min_length: int minimum length peptide to generate. Default=3
'''
def tryptic(sequences, number_digests, miss_prob=0, save_dir='./', save_name='digestion.tsv', min_length=3):
    save_dir = __make_valid_dir_string(save_dir)
    __make_dir(save_dir)

    to_digest = []
    if number_digests <= len(sequences):
        for _ in range(number_digests):
            to_digest.append(choice(sequences))
    else:
        to_digest += sequences
        for _ in range(number_digests - len(sequences)):
            to_digest.append(choice(sequences))

    peptides = []
    o = open(save_dir + save_name, 'w')
    form = '{}\t{}\t{}\t{}\n'
    o.write(form.format('peptide', 'parent-name', 'parent-sequence', 'start-location'))

    for digest in to_digest:
        seq = digest['sequence']
        name = digest['name']
        this_pep, start = __tryptic_digest(seq, miss_prob)
        #ensure that no peptide is shorter than the minimum length
        while len(this_pep) < min_length:
            this_pep, start = __tryptic_digest(seq, miss_prob)

        peptides.append(this_pep)
        o.write(form.format(this_pep, name, seq, start))

    return peptides
