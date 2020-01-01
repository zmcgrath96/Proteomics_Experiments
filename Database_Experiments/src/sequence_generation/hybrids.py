from random import randint, choice
from math import ceil
from utils import __make_valid_csv_file, __make_valid_dir_string, __file_exists

hybrid_peptide_file_name = 'hybrid_peptides.csv'
hybrid_protein_file_name = 'hybrid_proteins.csv'

'''__make_hybrid_pep

DESC:
    generate a hybrid peptide
PARMS:
    hybrid_prot: a dictionary with entries 
        {
            hybrid_protein: str,
            left_parent_end: int,
            name: str,
        }
OPTIONAL:
    min_length: int minimum length peptide to generate. Default=4
    max_length: int maximum length peptide to generate. Default=20
RETURNS:
    dictionary of form
    {
        hybrid_peptide_sequence: str,
        hybrid_parent_name: str,
        hybrid_parent_sequence: str,
        start_index: int, 
        end_index: int
    }
'''
def __make_hybrid_pep(hybrid_prot, min_length=4, max_length=20):
    l = randint(min_length, max_length)
    left_contr = randint(1, l-1)
    right_contr = l - left_contr
    j_site = hybrid_prot['left_parent_end']
    pep = hybrid_prot['protein'][j_site-left_contr:j_site+right_contr]
    return {
        'peptide_sequence': pep,
        'parent_name': hybrid_prot['name'],
        'parent_sequence': hybrid_prot['protein'],
        'start_index': hybrid_prot['protein'].index(pep),
        'end_index': hybrid_prot['protein'].index(pep) + len(pep)
    }


'''generate_peptides

DESC:
    generate hybrid peptides. Guaranteed to capture the junction
PARAMS:
    hybrid_prots: list of dictionaries. Should be of the form 
        [{
            left_parent_name: str,
            right_parent_name: str,
            left_parent_sequence: str,
            right_parent_sequence: str,
            left_parent_start: int,
            right_parent_start: int, 
            left_parent_contribution: int,
            right_parent_contribution: int,
            hybrid_protein: str, 
            name: str
        }]
OPTIONAL:
    num_gen: int number of hybrid peptides to genereate. if num_gen < len(hybrid_prots), then 
        the remaining will be ignored. After one is generated from. Default=10
    peptide_name_prefix: str prefix to add as the name to each protein. Default=HYBRID_PEPTIDE
    min_length: int minimum length peptide to generate. Default=4
    max_length: int maximum length peptide to generate. Default=20
RETURNS:
    l1, l2
    l1: list of dictionaries of the form 
        [{
            hybrid_peptide_sequence: str,
            hybrid_peptide_name: str,
            hybrid_parent_name: str,
            hybrid_parent_sequence: str,
            start_index: int,
            end_index: int
        }]
    l2: list of hybrid_proteins of the input form
'''
def generate_peptides(hybrid_prots, num_gen=10, peptide_name_prefix='HYBRID_PEPTIDE_', min_length=4, max_length=20):
    first_round = num_gen if num_gen <= len(hybrid_prots) else len(hybrid_prots)
    second_round = num_gen if num_gen > len(hybrid_prots) else 0
    hybrid_peps = []
    name_c = 0
    for i in range(first_round):
        print('Generating hybrid protein {}/{}[{}%]\r'.format(name_c, num_gen, int(float(name_c)/float(num_gen) * 100)), end="")
        hyb_pep = __make_hybrid_pep(hybrid_prots[i], min_length=min_length, max_length=max_length)
        hyb_pep['peptide_name'] = peptide_name_prefix + str(name_c).zfill(ceil(num_gen/10))
        hybrid_peps.append(hyb_pep)
        name_c += 1

    if second_round == 0: 
        print('Finshed generating hybrid peptides')
        return hybrid_peps, hybrid_prots

    for _ in range(second_round):
        print('Generating hybrid protein {}/{}[{}%]\r'.format(name_c, num_gen, int(float(name_c)/float(num_gen) * 100)), end="")
        hyb_pep = __make_hybrid_pep(choice(hybrid_prots), min_length=min_length, max_length=max_length)
        hyb_pep['peptide_name'] = peptide_name_prefix + str(name_c).zfill(ceil(num_gen/10))
        hybrid_peps.append(hyb_pep)
        name_c += 1

    print('Finished generating hybrid peptides')
    return hybrid_peps, hybrid_prots

'''generate_proteins

DESC:
    Generate a number of parent proteins with certain constraints
PARMS:
    num_gen: int number of hybrid_proteins to generate
    prots: list of dictionaries of the form
        [{
            'name': str,
            'sequence': str
        }]
OPTIONAL:
    min_contribution: int minimum number of AAs to use from each parent. Default=10
    name_prefix: str prefix to add before every name of protein. Default=HYBRID_
RETURNS:
    list of dictionaries. From is
    [{
        left_parent_name: str,
        right_parent_name: str,
        left_parent_sequence: str,
        right_parent_sequence: str,
        left_parent_end: int,
        right_parent_start: int, 
        left_parent_contribution: int,
        right_parent_contribution: int,
        hybrid_protein: str, 
        name: str, 
    }]
'''
def generate_proteins(num_gen, prots, min_contribution=10, name_prefix='HYBRID_'):
    hybrids = []
    for i in range(num_gen):
        print('Generating hybrid protein {}/{}[{}%]\r'.format(i, num_gen, int(float(i)/float(num_gen) * 100)), end="")
        left_parent = prots[randint(0, len(prots)-1)]
        right_parent = prots[randint(0, len(prots)-1)]
        left_end = randint(0, len(left_parent['sequence'])-1)
        right_start = randint(0, len(right_parent['sequence'])-1)
        while left_end < min_contribution:
            left_end = randint(0, len(left_parent['sequence'])-1)
        while len(right_parent['sequence']) - right_start < min_contribution:
            right_start = randint(0, len(right_parent['sequence'])-1)
        hybrid = left_parent['sequence'][:left_end] + right_parent['sequence'][right_start:]
        name = name_prefix + str(i).zfill(ceil(num_gen/10))
        hybrid_d = {
            'left_parent_name': left_parent['name'],
            'right_parent_name': right_parent['name'],
            'left_parent_sequence': left_parent['sequence'],
            'right_parent_sequence': right_parent['sequence'],
            'left_parent_end': left_end,
            'right_parent_start': right_start,
            'left_parent_contribution': left_end - 1,
            'right_parent_contribution': len(right_parent) - right_start,
            'protein': hybrid,
            'name': name
        }
        hybrids.append(hybrid_d)
    print('Finished generating hybrid proteins')
    return hybrids

'''save_peptides

DESC:
    save all the hybrid peptides to file in csv
PARAMS:
    hybrid_peps: list of dictionaries from the function generate peptides
OPTIONAL:
    save_dir: str path to the directory to save the file. Default=./
RETURNS:
    str the file name of the saved file
'''
def save_peptides(hybrid_peps, save_dir='./'):
    save_dir = __make_valid_dir_string(save_dir)
    f = save_dir + hybrid_peptide_file_name
    __make_valid_csv_file(f)
    form = '{},{},{},{},{}\n'
    with open(f, 'w') as o:
        o.write(form.format('peptide_sequence', 'peptide_name', 'parent_name', 'parent_sequence', 'start_index', 'end_index'))
        for pep in hybrid_peps:
            o.write(form.format(pep['peptide_sequence'], pep['peptide_name'], pep['parent_name'], pep['parent_sequence'], pep['start_index'], pep['end_index']))
    return f

'''save_proteins

DESC:
    save all the hybrid proteins to file in csv
PARAMS:
    hybrid_prots: list of proteins from the function generate proteins
OPTIONAL:
    save_dir: str path to the directory to save the file. Default=./
RETURNS:
    str the file name of the saved file
'''
def save_proteins(hybrid_prots, save_dir='./'):
    save_dir = __make_valid_dir_string(save_dir)
    f = save_dir + hybrid_protein_file_name
    __make_valid_csv_file(f)
    form = '{},{},{},{},{},{},{},{},{},{}\n'
    with open(f, 'w') as o:
        o.write(form.format('name', 'protein', 'left_parent_name', 'right_parent_name', 'left_parent_sequence', 'right_parent_sequence', 'left_parent_end', 'right_parent_start', 'left_parent_contribution', 'right_parent_contribution'))
        for p in hybrid_prots:
            o.write(form.format(p['name'], p['protein'], p['left_parent_name'], p['right_parent_name'], p['left_parent_sequence'], p['right_parent_sequence'], p['left_parent_end'], p['right_parent_start'], p['left_parent_contribution'], p['right_parent_contribution']))
    return f

'''read_peptides

DESC:
    read in a csv with hybrid peptides into memory
PARAMS:
    hybrid_pep_file: str path to a csv file with the hybrid peptides
RETURNS:
    list of dictionaries of the form
    [{
        hybrid_peptide_sequence: str,
        hybrid_peptide_name: str,
        hybrid_parent_name: str,
        hybrid_parent_sequence: str,
        start_index: int,
        end_index: int
    }]
'''
def read_peptides(hybrid_pep_file):
    peps = []
    read_header = False
    if not __file_exists(hybrid_pep_file):
        raise Exception('Could not find file {}'.format(hybrid_pep_file))
    with open(hybrid_pep_file, 'r') as i:
        for line in i:
            if not read_header:
                read_header = True
                continue
            l = line.split(',')
            p = {
                'peptide_sequence': l[0],
                'peptide_name': l[1],
                'parent_name': l[2],
                'parent_sequence': l[3],
                'start_index': l[4], 
                'end_index': l[5]
            }
            peps.append(p)
    return peps

'''read_proteins

DESC:
    read in a csv with hybrid proteins into memory
PARAMS:
    hybrid_prot_file: str path to a csv file with the hybrid proteins
RETURNS:
    list of dictionaries of the form
    [{
        left_parent_name: str,
        right_parent_name: str,
        left_parent_sequence: str,
        right_parent_sequence: str,
        left_parent_end: int,
        right_parent_start: int, 
        left_parent_contribution: int,
        right_parent_contribution: int,
        hybrid_protein: str, 
        name: str, 
    }]
'''
def read_proteins(hybrid_prot_file):
    prots = []
    read_header = False
    if not __file_exists(hybrid_prot_file):
        raise Exception('Could not find file {}'.format(hybrid_prot_file))
    with open(hybrid_prot_file, 'r') as i:
        for line in i:
            if not read_header:
                read_header = True
                continue
            l = line.split(',')
            p = {
                'left_parent_name': l[2],
                'right_parent_name': l[3],
                'left_parent_sequence': l[4],
                'right_parent_sequence': l[5],
                'left_parent_end': l[6],
                'right_parent_start': l[7], 
                'left_parent_contribution': l[8],
                'right_parent_contribution': l[9],
                'protein': l[1], 
                'name': l[0], 
            }
            prots.append(p)
    return prots