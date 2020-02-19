from sequences.digest import tryptic, random_digest
from math import ceil
from random import choice, randint

digest_functions = {
    'random': random_digest,
    'trypsin': tryptic
}

############################################################################################################
#           PRIVATE FUNCTIONS
############################################################################################################
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

'''__make_hybrid_peps_brute_force

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
    max_contribution: int max contribution from each side. Default=10
    min_length: mininum length peptide to create. Default=2
RETURNS:
    list of dictionarie of form
    [{
        hybrid_peptide_sequence: str,
        hybrid_parent_name: str,
        hybrid_parent_sequence: str,
        start_index: int, 
        end_index: int
    }]
'''
def __make_hybrid_peps_brute_force(hybrid_prot, max_contribution=10, min_length=2):
    hyb_peps = []
    j_site = hybrid_prot['left_parent_end']

    for i in range(1, max_contribution):
        for j in range(1, max_contribution):
            pep = hybrid_prot['protein'][j_site - i:j_site + j]
            if len(pep) < min_length:
                continue
            hyb_peps.append({
                'peptide_sequence': pep,
                'parent_name': hybrid_prot['name'],
                'parent_sequence': hybrid_prot['protein'],
                'start_index': hybrid_prot['protein'].index(pep),
                'end_index': hybrid_prot['protein'].index(pep) + len(pep),
            })
    return hyb_peps

'''__generate_hybrids

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
    list of dictionaries of the form 
        {
            hybrid_peptide_sequence: str,
            hybrid_peptide_name: str,
            hybrid_parent_name: str,
            hybrid_parent_sequence: str,
            start_index: int,
            end_index: int
        }
'''
def __generate_hybrids(hybrid_prots, num_gen=10, peptide_name_prefix='HYBRID_PEPTIDE_', min_length=4, max_length=20):
    first_round = num_gen if num_gen <= len(hybrid_prots) else len(hybrid_prots)
    second_round = num_gen if num_gen > len(hybrid_prots) else 0
    hybrid_peps = []
    name_c = 0

    fill_zeros = len(str(ceil(num_gen / 10)))

    for i in range(first_round):
        print('Generating hybrid peptide {}/{}[{}%]\r'.format(name_c, num_gen, int(float(name_c)/float(num_gen) * 100)), end="")
        hyb_pep = __make_hybrid_pep(hybrid_prots[i], min_length=min_length, max_length=max_length)
        hyb_pep['peptide_name'] = peptide_name_prefix + str(name_c).zfill(fill_zeros)
        hybrid_peps.append(hyb_pep)
        name_c += 1

    if second_round == 0: 
        print('\nFinshed generating hybrid peptides')
        return hybrid_peps

    for _ in range(second_round):
        print('Generating hybrid peptide {}/{}[{}%]\r'.format(name_c, num_gen, int(float(name_c)/float(num_gen) * 100)), end="")
        hyb_pep = __make_hybrid_pep(choice(hybrid_prots), min_length=min_length, max_length=max_length)
        hyb_pep['peptide_name'] = peptide_name_prefix + str(name_c).zfill(fill_zeros)
        hybrid_peps.append(hyb_pep)
        name_c += 1

    print('\nFinished generating hybrid peptides')
    return hybrid_peps

'''__generate_hybrids_brute_force

DESC:
    Generate all possible hybrid peptides from a window of a hybrid 
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
    peptide_name_prefix: str prefix to add as the name to each protein. Default=HYBRID_PEPTIDE
    min_length: int minimum length peptide to generate. Default=2
    max_contribution: int maximum to allow each side of the hybrid to allow. Default=10
RETURNS:
    list of dictionaries of the form 
        {
            hybrid_peptide_sequence: str,
            hybrid_peptide_name: str,
            hybrid_parent_name: str,
            hybrid_parent_sequence: str,
            start_index: int,
            end_index: int
        }
'''
def __generate_hybrids_brute_force(hybrid_prots, peptide_name_prefix='HYBRID_PEPTIDE', min_length=2, max_contribution=10):
    hybrid_peps = []
    name_c = 0

    for prot in hybrid_prots:
        hybrid_peps += __make_hybrid_peps_brute_force(prot, max_contribution=max_contribution, min_length=min_length)
    fill_zeros = len(str(ceil(len(hybrid_peps) / 10)))
    for pep in hybrid_peps:
        pep['peptide_name'] = peptide_name_prefix + str(name_c).zfill(fill_zeros)
        name_c += 1
    
    return hybrid_peps

############################################################################################################
#          END PRIVATE FUNCTIONS
############################################################################################################

'''gen_peptides

DESC:
    generates peptides from proteins given
PARAMS:
    proteins: list of dictionaries of the form {'name': str, 'sequence': str}
              NOTE: if hybrid_list set to true, form is 
              {
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
              }
    n: number of peptides to generate
OPTIONAL:
    min_length: int minimum length peptide to generate. Default=3
    max_length: int maximum length peptide to generate. Default=20
    digest: str type of digest to perform. Default=random
    hybrid_list: bool if the proteins passed in are hybrids and you wish to capture a junction point, set to True. Default=False
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
def gen_peptides(proteins, n, min_length=3, max_length=20, digest='random', hybrid_list=False):
    digest = 'random' if digest.lower() not in digest_functions.keys() else digest.lower()
    
    if not hybrid_list:
        return digest_functions[digest](proteins, n, min_length=min_length, max_length=max_length)

    else:
        # return __generate_hybrids(proteins, num_gen=n, min_length=min_length, max_length=max_length)
        return __generate_hybrids_brute_force(proteins, min_length=min_length)


    