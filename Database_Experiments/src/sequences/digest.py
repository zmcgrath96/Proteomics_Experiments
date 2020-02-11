from random import randint, choice, random, shuffle
from utils import __make_dir, __make_valid_dir_string
from math import ceil
from copy import deepcopy
####################################################################
#                   PRIVATE FUNCTIONS
####################################################################
'''__verify_length

DESC:
    make sure proteins are long enough
PARAMS: 
    seqs: list of dictionaries {'name': str, 'sequence': str}
    min_length: int 
RETURNS:
    list of dictionaries of the same form 
'''
def __verify_length(seqs, min_length):
    verified = []
    for seq in seqs:
        if len(seq['sequence']) >= min_length:
            verified.append(seq)
    return verified
####################################################################
#                  END PRIVATE FUNCTIONS
####################################################################

####################################################################
#                   TRYPTIC DIGEST FUNCTIONS
####################################################################
'''__tryptic_digest

DESC:
    Perfrom the actual digestion 
PARAMS: 
    sequence: full length sequence to perform digestion on
OPTIONAL:
    miss_prob: 
    rand_cut_prob: float value [0, 1] probability of cutting peptide short randomly. Default=.05
'''
def __tryptic_digest(sequence, miss_prob=0, rand_cut_prob=.05):
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
    
    for i in range(len(this_pep)-2, 0, -1):
        if random() < rand_cut_prob:
            this_pep = this_pep[i:]

    return this_pep, sequence.index(this_pep)

'''tryptic

DESC:
    Generate peptides with tryptic digest from sequences with missed cleavages at a probability
PARAMS:
    sequences: list of dictionary objects with entries {'sequence': string, 'name': string}
    number_digests: int number of peptides to generate
OPTIONAL:
    miss_prob: float the probability of a missed cleavage. Should range [0, 1). Default=0
    save_dir: string the directory in which to save the digestion file. Default=./
    save_name: string the name to save the digestion information in. Default=peptides.tsv
    min_length: int minimum length peptide to generate. Default=3
    rand_cut_prob: float value [0, 1] probability to randomly cut a peptide short. Default=.05
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
def tryptic(sequences, number_digests, peptide_prefix='peptide_', miss_prob=0, save_dir='./', save_name='peptides.tsv', min_length=3, max_length=20, rand_cut_prob=.05):
    save_dir = __make_valid_dir_string(save_dir)
    __make_dir(save_dir)

    sequences = __verify_length(sequences, min_length)

    to_digest = []
    if number_digests <= len(sequences):
        for _ in range(number_digests):
            to_digest.append(choice(sequences))
    else:
        to_digest += sequences
        for _ in range(number_digests - len(sequences)):
            to_digest.append(choice(sequences))

    fill_zeros = len(str(ceil(number_digests / 10)))

    peptides = []
    o = open(save_dir + save_name, 'w')
    form = '{}\t{}\t{}\t{}\t{}\t{}\n'
    o.write(form.format('peptide_name', 'peptide', 'parent_name', 'parent_sequence', 'start_location', 'end_location'))

    digest_count = 0
    for digest in to_digest:
        print('Creating peptide from digest {}/{}[{}%]\r'.format(digest_count, number_digests, int(float(digest_count)/float(number_digests) * 100)), end="")
        seq = digest['sequence']
        name = digest['name']

        # make sure the protein is long enough to operate on
        if len(seq) < min_length: 
            print('\nProtein sequence too short. Name: {}\t sequence: {}'.format(name, seq))
            continue

        pep_name = peptide_prefix + str(digest_count).zfill(fill_zeros)
        this_pep, start = __tryptic_digest(seq, miss_prob=miss_prob, rand_cut_prob=rand_cut_prob)
        #ensure that no peptide is shorter than the minimum length
        if len(this_pep) < min_length:
            idx = seq.index(this_pep)
            if idx - (min_length - len(this_pep)) >= 0:
                this_pep = seq[idx - (min_length - len(this_pep)): idx + min_length]
                start = idx
            else:
                this_pep = seq[idx: idx + min_length]
                start = idx

        end = start + len(this_pep)
        # if the peptide is too long, cut from the left side
        if len(this_pep) > max_length:
            start_pep = len(this_pep) - max_length
            this_pep = this_pep[start_pep:]
            start = seq.index(this_pep)

        pep_obj = {
            'peptide_name': pep_name,
            'peptide_sequence': this_pep,
            'parent_name': name,
            'parent_sequence': seq,
            'start_index': start, 
            'end_index': end
        }
        peptides.append(pep_obj)
        o.write(form.format(pep_name, this_pep, name, seq, start, end))
        digest_count += 1

    print('\nFinished digestion')
    return peptides

####################################################################
#               END TRYPTIC DIGEST FUNCTIONS
####################################################################

####################################################################
#                     RANDOM DIGEST FUNCTIONS
####################################################################

'''random_digest

DESC:
PARAMS:
    proteins: list of dictionaries of the form {'sequence': string, 'name': string}
    n: number of peptides to create
OPTIONAL:
    peptide_prefix: str prefix to give all peptide names. Default='peptide_'
    min_length: int minimum length peptide to create. Default=3
    max_length: int maximum length peptide to create. Default=20
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
def random_digest(proteins, n, peptide_prefix='peptide_', min_length=3, max_length=20):
    digested = []
    fill_zeros = len(str(ceil(n / 10)))

    to_digest = deepcopy(proteins)
    if n > len(proteins):
        for _ in range(len(proteins), n):
            to_digest.append(proteins[randint(0, len(proteins) - 1)])

    for i in range(min(n, len(to_digest))):
        seq = to_digest[i]['sequence']
        r = randint(min_length, max_length)
        start = randint(0, len(seq) - max_length)
        pep = seq[start : start + r]
        pep_name = peptide_prefix + str(i).zfill(fill_zeros)

        d = {
            'peptide_name': pep_name,
            'peptide_sequence': pep,
            'parent_name': to_digest[i]['name'],
            'parent_sequence': seq,
            'start_index': start,
            'end_index': start + r
        }
        digested.append(d)

    return digested

####################################################################
#                   END RANDOM DIGEST FUNCTIONS
####################################################################