from random import randint, choice, random
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
        'start_index': int, 
        'end_index': int
    }
'''
def load_digest(digest_file):
    digests = {}
    with open(digest_file, 'r') as o:
        for i, line in enumerate(o):
            if i == 0: # skip the header line
                continue
            l = line.split('\t')
            num = '0' + str(i-1) if i < 11 else str(i-1)
            name = 'peptide_' + num
            entry = {'peptide_name': l[0], 'peptide_sequence': l[1], 'parent_name': l[2], 'parent_sequence': l[3], 'start_index': int(l[4]), 'end_index': int(l[5])}
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
    
    return this_pep, sequence.index(this_pep)

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
def tryptic(sequences, number_digests, peptide_prefix='peptide_', miss_prob=0, save_dir='./', save_name='digestion.tsv', min_length=3, max_length=20):
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
    form = '{}\t{}\t{}\t{}\t{}\t{}\n'
    o.write(form.format('peptide-name', 'peptide', 'parent-name', 'parent-sequence', 'start-location', 'end-location'))

    digest_count = 0
    for digest in to_digest:
        seq = digest['sequence']
        name = digest['name']
        pep_name = peptide_prefix + str(digest_count) if digest_count > 9 else peptide_prefix + '0' + str(digest_count)
        this_pep, start = __tryptic_digest(seq, miss_prob)
        end = start + len(this_pep)
        #ensure that no peptide is shorter than the minimum length
        while len(this_pep) < min_length:
            this_pep, start = __tryptic_digest(seq, miss_prob)
        
        # if the peptide is too long, cut from the left side
        if len(this_pep) > max_length:
            start_pep = len(this_pep) - max_length
            this_pep = this_pep[start_pep:]
            start = seq.index(this_pep)

        peptides.append(this_pep)
        o.write(form.format(pep_name, this_pep, name, seq, start, end))
        digest_count += 1

    return peptides
