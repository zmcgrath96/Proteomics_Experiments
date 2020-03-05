from numbers import Number
from utils import __make_valid_dir_string, __make_dir
from file_io import JSON

######################################################################
#                   CONSTANTS
######################################################################

from numbers import Number

DIVIDER_WIDTH = 70
SUMMARY_FILE_NAME = 'summary.txt'
SUMMARY_JSON_FILE_NAME = 'summary.json'
SECTION_DIVIDER = '\n\n' + '=' * DIVIDER_WIDTH + '\n\n\n\n\n'
HEADER_UNDERLINE = '\n' + '-' * DIVIDER_WIDTH + '\n'
PROTEIN_SUMMARY_HEADER = 'EXPERIMENT PROTEIN INFORMATION'
PEPTIDE_SUMMARY_HEADER = 'EXPERIMENT PEPTIDE INFORMATION'
NON_HYBRID_SUMMARY_HEADER = 'NON HYBRID PEPTIDE SUMMARY'
HYBRID_SUMMARY_HEADER = 'HYBRID PEPTIDE SUMMARY'
LEFT_PARENT_SUMMARY_HEADER = 'LEFT PARENT SUMMARY'
RIGHT_PARENT_SUMMARY_HEADER = 'RIGHT PARENT SUMMARY'
HEADER_ROW_NAMES_PROTEIN_SUMMARY = [
    'number of proteins',           # int
    'number non-hybrid proteins',   # int
    'number hybrid proteins'        # int
]
HEADER_ROW_NAMES_PEPTIDE_SUMMARY = [
    'number of peptides',           # int
    'number non-hybrid peptides',   # int
    'number hybrid peptides'        # int
]
HEADER_ROW_NAMES_EXPERIMENT_NON_HYBRID_SUMMARY = [
    'total peptides ',               # int
    'correct predictions ',          # int
    '(%) ',                          # float
    'correct parent protein ',       # int
    '(%) ',                          # float
    'near miss predictions* ',       # float
    '(%) ',                          # float
    'correct starting position* ',   # int
    '(%) ',                          # float
    'correct peptide length* ',      # int
    '(%) '                           # float
]
HEADER_ROW_NAMES_EXPERIMENT_HYBRID_SUMMARY = [
    'total peptides ',               # int
    'correct predictions** ',        # int
    '(%) ',                          # float
    'correct predictions*** ',       # int
    '(%) ',                          # float                              
    'near miss**** ',                # int
    '(%) ',                          # float
    'correct protein ',              # int
    '(%) ',                          # float
    'correct starting position ',    # int
    '(%) ',                          # float
    'correct length ',               # int
    '(%) '                           # float
]

__cols = lambda x: ''.join(['{}\t' for _ in range(x)])

HEADER_ROW_PROTEIN_SUMMARY = __cols(len(HEADER_ROW_NAMES_PROTEIN_SUMMARY)).format(*HEADER_ROW_NAMES_PROTEIN_SUMMARY)
HEADER_ROW_PEPTIDE_SUMMARY = __cols(len(HEADER_ROW_NAMES_PEPTIDE_SUMMARY)).format(*HEADER_ROW_NAMES_PEPTIDE_SUMMARY)
HEADER_ROW_EXPERIMENT_NON_HYBRID_SUMMARY = __cols(len(HEADER_ROW_NAMES_EXPERIMENT_NON_HYBRID_SUMMARY)).format(*HEADER_ROW_NAMES_EXPERIMENT_NON_HYBRID_SUMMARY)

NEAR_MISS_CONTEXT = '\n\n* near miss \n \
A near miss means that the predicted peptide sequence identified the correct \n \
parent protein and was ONLY ONE of the following:\n \
    (1)  Predicted length was incorrect\n \
    (2)  The starting position was incorrect\n'

HYBRID_CORRECT_PREDICTIONS_CONTEXT = '\n\n** correct hybrid prediction\n \
A correct hybrid prediction is very unlikely. For it to be completely correct, \n \
it must have done the following: \n \
    (1)  Had predicted both parent proteins in the top 2 rankings \n \
    (2)  For both of these parents, it must have \n \
        (a)  Correctly identified the starting position \n \
        (b)  Correctly identified the length of the protein \n'

HYBRID_CORRECT_PARENT_PREDICTIONS_CONTEXT = '\n\n*** hybrid correct parent prediction \n \
This is similar to the non-hybrid parent prediction. A parent protein of the hybrid \n \
was correctly identified. This means that \n \
    (1)  The correct parent protein was reported \n \
    (2)  The correct starting position within the protein was reported \n \
    (3)  The correct length of the parent was predicted'

HYBRID_NEAR_MISS_CONTEXT = '\n\n* hybrid near miss \n \
A hybrid near miss means both parent proteins were correctly identified\n \
and each parent protein was ONLY ONE of the following: \
    (1)  Predicted length was incorrect\n \
    (2)  The starting position was incorrect\n'

def __pad_and_center(num: Number, s: str) -> str:
    '''
    pad and center a number in a string
    
    Inputs:
        num: the number to center
        s:   the string to center the number int
    Outputs:
        centered: string representing the number
    '''
    if len(str(num)) >= len(s):
        return str(num) 
    
    center = ''
    side_len = (len(s) - len(str(num))) // 2 
    center += ' ' * side_len + str(num) + ' ' * side_len
    return center


######################################################################
#                         END CONSTANTS
######################################################################

######################################################################
#                       PRIVATE FUNCTIONS
######################################################################

def __find_pep(peptide_name: str, exp: dict) -> dict:
    '''
    Find the peptide information in the experiment header
    
    Inputs:
        peptide_name:    string name of the peptide
        exp:             dictionary of the experiment data
    Outputs:
        dictionary corresponding to the correct peptide. {} if not found
    '''
    for pep in exp['experiment_info']['peptides']:
        if peptide_name.lower() == pep['peptide_name'].lower():
            return pep
    return {}

def __find_prot(protein_name: str, exp: dict) -> dict:
    '''
    Find the protein information in the experiment header
    
    Inputs:
        protein_name:    string name of the protein
        exp:             dictionary of the experiment data
    Outputs:
        dictionary corresponding to the correct protein. {} if not found
    '''
    for prot in exp['experiment_info']['proteins']:
        if protein_name.lower() == prot['name'].lower():
            return prot
    return {}

def __prediction_stats_non_hyb(exp: dict) -> tuple:
    '''
    Find prediction stats of correct predictions, near misses, correct parents, correct starting position, correct length
    
    Inputs:
        exp:    dictionary of the experiment data
    Outputs:
        tuple of stats in the form (totally correct, near miss, correct parents, correct starting position, correct length)
    '''
    correct_count = 0
    near_miss = 0
    correct_parents = 0
    correct_starting_position = 0
    correct_length = 0
    
    for pep_name, pep in exp['experiment'].items():
        if 'hybrid' in pep_name.lower():
            continue

        full_correct = False
        # find the correct length, sequence, parent protein
        peptide_info = __find_pep(pep_name, exp)              # known information
        peptide_prediction = pep['analysis']['sequence_predictions'][0]  # prediction information

        corr_par = peptide_info['parent_name'] == peptide_prediction['protein_name']
        corr_start_pos = peptide_info['start_index'] == peptide_prediction['starting_pos']
        corr_len = peptide_info['end_index'] - peptide_info['start_index'] == peptide_prediction['predicted_length']

        if corr_par and corr_start_pos and corr_len:
            full_correct = True
            correct_count += 1

        if corr_par and (corr_start_pos or corr_len) and not full_correct:
            near_miss += 1

        if corr_par and corr_start_pos and not corr_len:
            correct_starting_position += 1

        if corr_par and not corr_start_pos and corr_len:
            correct_length += 1

        if corr_par:
            correct_parents += 1
    
    return (correct_count, near_miss, correct_parents, correct_starting_position, correct_length)

def __hybrid_single_parent_prediction(sub_peptide_info: dict, candidates: list) -> tuple:
    '''
    Try and match one of the potential candidates to one of the known parents
    
    Inputs:
        sub_peptide_info:   Information on the half of the peptide we're trying to match
                        {
                            'parent_name': str,
                            'starting_position': int,
                            'length': int
                        }
        candidates:          List of potential candidates (the top scoring proteins).
                             Ideally limited to 2 to replicate "real" data
    Outputs:
        tuple of the following booleans
            (fully correct, correct parent, near miss, correct starting position, correct length)
    '''
    fully_corr, corr_parent, near_miss, corr_start_pos, corr_len = False, False, False, False, False
    
    # determine if there is a potential match at all
    candidates = [x for x in candidates if x['protein_name'] == sub_peptide_info['parent_name']]
    if len(candidates) == 0:
        return (False, False, False, False, False)
    corr_parent = True
    
    # determine if there is a starting position that matches
    corr_start_pos_candidates = [x for x in candidates if x['starting_pos'] == sub_peptide_info['starting_position']]
    if len(corr_start_pos_candidates) > 0:
        corr_start_pos = True
        candidates = corr_start_pos_candidates
        
    # determine if there is a correct_len
    corr_len_candidates = [x for x in candidates if x['predicted_length'] == sub_peptide_info['length']]
    if len(corr_len_candidates) > 0:
        corr_len = True
        
    fully_corr = corr_parent and corr_start_pos and corr_len
    near_miss = corr_parent and (corr_start_pos or corr_len) and not fully_corr
    
    return (fully_corr, corr_parent, near_miss, corr_start_pos, corr_len)
    
def __get_sub_peptide_info(peptide_name: str, exp: dict) -> (dict, dict):
    '''
    Get the left and right sub-peptide info from the hybrid
    
    Inputs:
        peptide_name:   string name of the peptide to split up
        exp:            dictionary of the experiment data
    Outputs:
        left_sub_peptide, rigth_sub_peptide
        both of these are dictionaries of 
        {
            'parent_name': str,
            'starting_position': int,
            'length': int
        }
    '''
    left_sub_peptide = {}
    right_sub_peptide = {}
    
    hyb_pep_info = __find_pep(peptide_name, exp)
    hyb_prot_info = __find_prot(hyb_pep_info['parent_name'], exp)
    
    left_sub_peptide['parent_name'] = hyb_prot_info['left_parent_name']
    right_sub_peptide['parent_name'] = hyb_prot_info['right_parent_name']
    
    left_sub_peptide['starting_position'] = hyb_pep_info['start_index'] # since left peptide indexed by hybrid and left starts, we're set here
    right_sub_peptide['starting_position'] = hyb_prot_info['right_parent_start'] # we know the junction will always be here
    
    left_sub_peptide['length'] = hyb_prot_info['left_parent_contribution'] - hyb_pep_info['start_index']
    right_sub_peptide['length'] = len(hyb_pep_info['peptide_sequence']) - left_sub_peptide['length']
    
    return left_sub_peptide, right_sub_peptide

def __prediction_stats_hyb(exp: dict) -> (int, dict, dict):
    '''
    Find prediction stats of correct predictions, near misses, correct parents, correct starting position, correct length
    
    Inputs:
        exp:    dictionary of the experiment data
    Outputs:
        correct_count, left_stats, right_stats
            correct_count:     integer number of totally correct hits
            left_stats:        dictionary of stats for the left parent
            right_stats:       dictionary of stats for the right parent
                stats are in form: {
                    'correct_count': int,
                    'near_miss': int,
                    'correct_parents': int,
                    'correct_starting_position': int,
                    'correct_length': int
                }
    '''
    keys = ['correct_count', 'near_miss', 'correct_parents', 'correct_starting_position', 'correct_length']
    left_parent = {}
    right_parent = {}
    for k in keys:
        left_parent[k] = 0
        right_parent[k] = 0
    correct_count = 0
    
    for pep_name, pep in exp['experiment'].items():
        if 'hybrid' not in pep_name.lower():
            continue
        candidates = pep['analysis']['sequence_predictions'][:2]
        # get each of the sub_peptide info
        left_sub_pep, right_sub_pep = __get_sub_peptide_info(pep_name, exp)
        left_stats = __hybrid_single_parent_prediction(left_sub_pep, candidates)
        right_stats = __hybrid_single_parent_prediction(right_sub_pep, candidates)
        
        # (fully_corr, corr_parent, near_miss, corr_start_pos, corr_len)
        left_parent['correct_count'] += 1 if left_stats[0] else 0
        right_parent['correct_count'] += 1 if right_stats[0] else 0
        left_parent['correct_parents'] += 1 if left_stats[1] else 0
        right_parent['correct_parents'] += 1 if right_stats[1] else 0
        left_parent['near_miss'] += 1 if left_stats[2] else 0
        right_parent['near_miss'] += 1 if right_stats[2] else 0
        left_parent['correct_starting_position'] += 1 if left_stats[3] else 0
        right_parent['correct_starting_position'] += 1 if right_stats[3] else 0
        left_parent['correct_length'] += 1 if left_stats[4] else 0
        right_parent['correct_length'] += 1 if right_stats[4] else 0
        
        correct_count += 1 if left_stats[0] and right_stats[0] else 0
        
    return correct_count, left_parent, right_parent
    

def __summary_header(prot_info: dict, pep_info: dict, summary_dict: dict) -> (str, dict):
    '''
    Create a string to return that has the basic experiment information in it
    
    Inputs:
        prot_info: dictionary of protein information
        pep_info:  dictionary of peptide information
        summary_dict: dictionary to save summary information in
    Outputs:
        header, summary_dict    
            header:    string of the header information
            summary_dict: same dictionary passed in with new information in it
    '''
    header = ''
    header += PROTEIN_SUMMARY_HEADER + HEADER_UNDERLINE + HEADER_ROW_PROTEIN_SUMMARY + '\n'
    
    prot_cols = __cols(len(HEADER_ROW_NAMES_PROTEIN_SUMMARY))
    pep_cols = __cols(len(HEADER_ROW_NAMES_PEPTIDE_SUMMARY))
    
    # protein summary is #, # non hybs, # hybs
    prot_count = len(prot_info)
    non_hyb_prot_count = len([x for x in prot_info if 'hybrid' not in x['name'].lower()])
    hyb_prot_count = prot_count - non_hyb_prot_count
    
    prot_count_s = __pad_and_center(prot_count, HEADER_ROW_NAMES_PROTEIN_SUMMARY[0])
    non_hyb_prot_count_s = __pad_and_center(non_hyb_prot_count, HEADER_ROW_NAMES_PROTEIN_SUMMARY[1])
    hyb_prot_count_s = __pad_and_center(hyb_prot_count, HEADER_ROW_NAMES_PROTEIN_SUMMARY[2])
    header += prot_cols.format(prot_count_s, non_hyb_prot_count_s, hyb_prot_count_s) + SECTION_DIVIDER + \
              PEPTIDE_SUMMARY_HEADER + HEADER_UNDERLINE + HEADER_ROW_PEPTIDE_SUMMARY +'\n'
    
    # peptide summary is #, # non hybs, # hybs
    pep_count = len(pep_info)
    non_hyb_pep_count = len([x for x in pep_info if 'hybrid' not in x['peptide_name'].lower()])
    hyb_pep_count = pep_count - non_hyb_pep_count
    
    pep_count_s = __pad_and_center(pep_count, HEADER_ROW_NAMES_PEPTIDE_SUMMARY[0])
    non_hyb_pep_count_s = __pad_and_center(non_hyb_pep_count, HEADER_ROW_NAMES_PEPTIDE_SUMMARY[1])
    hyb_pep_count_s = __pad_and_center(hyb_pep_count, HEADER_ROW_NAMES_PEPTIDE_SUMMARY[2])
    header += pep_cols.format(pep_count_s, non_hyb_pep_count_s, hyb_pep_count_s) + SECTION_DIVIDER
    
    # add info to the summary dict
    summary_dict['header'] = {
        'protein_info': {
            'count': prot_count,
            'non-hybrid_count': non_hyb_prot_count,
            'hybrid_count': hyb_prot_count
        },
        'peptide_info': {
            'count': pep_count, 
            'non-hybrid_count': non_hyb_pep_count,
            'hybrid_count': hyb_pep_count
        }
    }
    
    return header, summary_dict

def __non_hybrid_peptide_summary(exp: dict, summary_dict: dict) -> (str, dict):
    '''
    Generate a summary of the predictions done for non-hybrid peptides
    
    Inputs:
        exp:     dictionary with the experiment data
        summary_dict: dictionary to save summary information in
    Outputs:
        summary, summary_dict 
            summary:    string of the summary
            summary_dict: same dictionary passed in with new information in it
    '''
    #total peptides, correct predictions, (%), near miss predictions (%), correct parent protein, (%), correct starting position, (%), correct peptide length, (%)
    summary = ''
    summary += NON_HYBRID_SUMMARY_HEADER + HEADER_UNDERLINE
    
    # do the number crunching
    non_hybrid_peptide_count = len([x for x in exp['experiment_info']['peptides'] if 'hybrid' not in x['peptide_name'].lower()])
    percent_of = lambda x, c: round((float(x)/float(c)) * 100, 2)
    
    corr_c, near_miss, corr_par, corr_start_pos, corr_len = __prediction_stats_non_hyb(exp)
    corr_c_p = percent_of(corr_c, non_hybrid_peptide_count)
    near_miss_p = percent_of(near_miss, non_hybrid_peptide_count)
    corr_par_p = percent_of(corr_par, non_hybrid_peptide_count)
    corr_start_pos_p = percent_of(corr_start_pos, non_hybrid_peptide_count)
    corr_len_p = percent_of(corr_len, non_hybrid_peptide_count)
    
    rjust_width = lambda s: DIVIDER_WIDTH - len(s)
    
    summary += HEADER_ROW_NAMES_EXPERIMENT_NON_HYBRID_SUMMARY[0] + ' '.rjust(rjust_width(HEADER_ROW_NAMES_EXPERIMENT_NON_HYBRID_SUMMARY[0]), '~') + str(non_hybrid_peptide_count) + '\n' +\
               HEADER_ROW_NAMES_EXPERIMENT_NON_HYBRID_SUMMARY[1] + ' '.rjust(rjust_width(HEADER_ROW_NAMES_EXPERIMENT_NON_HYBRID_SUMMARY[1]), '~') + str(corr_c) + '\n' +\
               HEADER_ROW_NAMES_EXPERIMENT_NON_HYBRID_SUMMARY[2] + ' '.rjust(rjust_width(HEADER_ROW_NAMES_EXPERIMENT_NON_HYBRID_SUMMARY[2]), '~') + str(corr_c_p) + '\n' +\
               HEADER_ROW_NAMES_EXPERIMENT_NON_HYBRID_SUMMARY[3] + ' '.rjust(rjust_width(HEADER_ROW_NAMES_EXPERIMENT_NON_HYBRID_SUMMARY[3]), '~') + str(corr_par) + '\n' +\
               HEADER_ROW_NAMES_EXPERIMENT_NON_HYBRID_SUMMARY[4] + ' '.rjust(rjust_width(HEADER_ROW_NAMES_EXPERIMENT_NON_HYBRID_SUMMARY[4]), '~') + str(corr_par_p) + '\n' +\
               HEADER_ROW_NAMES_EXPERIMENT_NON_HYBRID_SUMMARY[5] + ' '.rjust(rjust_width(HEADER_ROW_NAMES_EXPERIMENT_NON_HYBRID_SUMMARY[5]), '~') + str(near_miss) +'\n' +\
               HEADER_ROW_NAMES_EXPERIMENT_NON_HYBRID_SUMMARY[6] + ' '.rjust(rjust_width(HEADER_ROW_NAMES_EXPERIMENT_NON_HYBRID_SUMMARY[6]), '~') + str(near_miss_p) + '\n' +\
               HEADER_ROW_NAMES_EXPERIMENT_NON_HYBRID_SUMMARY[7] + ' '.rjust(rjust_width(HEADER_ROW_NAMES_EXPERIMENT_NON_HYBRID_SUMMARY[7]), '~') + str(corr_start_pos) + '\n' +\
               HEADER_ROW_NAMES_EXPERIMENT_NON_HYBRID_SUMMARY[8] + ' '.rjust(rjust_width(HEADER_ROW_NAMES_EXPERIMENT_NON_HYBRID_SUMMARY[8]), '~') + str(corr_start_pos_p) + '\n' +\
               HEADER_ROW_NAMES_EXPERIMENT_NON_HYBRID_SUMMARY[9] + ' '.rjust(rjust_width(HEADER_ROW_NAMES_EXPERIMENT_NON_HYBRID_SUMMARY[9]), '~') + str(corr_len) + '\n' +\
               HEADER_ROW_NAMES_EXPERIMENT_NON_HYBRID_SUMMARY[10] + ' '.rjust(rjust_width(HEADER_ROW_NAMES_EXPERIMENT_NON_HYBRID_SUMMARY[10]), '~') + str(corr_len_p)
    
    summary += NEAR_MISS_CONTEXT + SECTION_DIVIDER
    
    # add info to the summary dict
    summary_dict['non-hybrid_peptide_summary'] = {
        'count': non_hybrid_peptide_count,
        'totally_correct': corr_c, 
        'totally_correct_percent': corr_c_p,
        'correct_parent': corr_par, 
        'correct_parent_percentage': corr_par_p,
        'near_miss': near_miss,
        'near_miss_percentage': near_miss_p,
        'correct_starting_position': corr_start_pos,
        'correct_starting_position_percentage': corr_start_pos_p,
        'correct_length': corr_len,
        'correct_length_percentage': corr_len_p
    }
    
    return summary, summary_dict
    
def __hybrid_peptide_summary(exp: dict, summary_dict: dict) -> (str, dict):
    '''
    Generate a summary of the predictions done for hybrid peptides
    
    Inputs:
        exp:     dictionary with the experiment data
        summary_dict: dictionary to save summary information in
    Outputs:
        summary, summary_dict 
            summary:    string of the summary
            summary_dict: same dictionary passed in with new information in it
    '''
    #total peptides, correct predictions, (%), near miss predictions (%), correct parent protein, (%), correct starting position, (%), correct peptide length, (%)
    summary = ''
    summary += HYBRID_SUMMARY_HEADER + HEADER_UNDERLINE
    
    hyb_pep_count = len([x for x in exp['experiment_info']['peptides'] if 'hybrid' in x['peptide_name'].lower()])
    
    percent_of = lambda x, c: round((float(x)/float(c)) * 100, 2)
    rjust_width = lambda s: DIVIDER_WIDTH - len(s)
    
    total_correct_count, left_info, right_info = __prediction_stats_hyb(exp)
    total_correct_count_p = percent_of(total_correct_count, hyb_pep_count)
    
    summary += HEADER_ROW_NAMES_EXPERIMENT_HYBRID_SUMMARY[0] + ' '.rjust(rjust_width(HEADER_ROW_NAMES_EXPERIMENT_HYBRID_SUMMARY[0]), '~') + str(hyb_pep_count) + '\n' + \
               HEADER_ROW_NAMES_EXPERIMENT_HYBRID_SUMMARY[1] + ' '.rjust(rjust_width(HEADER_ROW_NAMES_EXPERIMENT_HYBRID_SUMMARY[1]), '~') + str(total_correct_count) + '\n' + \
               HEADER_ROW_NAMES_EXPERIMENT_HYBRID_SUMMARY[2] + ' '.rjust(rjust_width(HEADER_ROW_NAMES_EXPERIMENT_HYBRID_SUMMARY[2]), '~') + str(total_correct_count_p) + '\n' 
    
    
    summary += '\n' + LEFT_PARENT_SUMMARY_HEADER + '\n'
    # correct predictions***', (%), near miss****, (%), correct protein, (%), correct starting position, (%), correct length, (%)
    summary += HEADER_ROW_NAMES_EXPERIMENT_HYBRID_SUMMARY[3] + ' '.rjust(rjust_width(HEADER_ROW_NAMES_EXPERIMENT_HYBRID_SUMMARY[3]), '~') + str(left_info['correct_count']) + '\n' + \
               HEADER_ROW_NAMES_EXPERIMENT_HYBRID_SUMMARY[4] + ' '.rjust(rjust_width(HEADER_ROW_NAMES_EXPERIMENT_HYBRID_SUMMARY[4]), '~') + str(percent_of(left_info['correct_count'], hyb_pep_count)) + '\n' + \
               HEADER_ROW_NAMES_EXPERIMENT_HYBRID_SUMMARY[5] + ' '.rjust(rjust_width(HEADER_ROW_NAMES_EXPERIMENT_HYBRID_SUMMARY[5]), '~') + str(left_info['near_miss']) + '\n' + \
               HEADER_ROW_NAMES_EXPERIMENT_HYBRID_SUMMARY[6] + ' '.rjust(rjust_width(HEADER_ROW_NAMES_EXPERIMENT_HYBRID_SUMMARY[6]), '~') + str(percent_of(left_info['near_miss'], hyb_pep_count)) + '\n' + \
               HEADER_ROW_NAMES_EXPERIMENT_HYBRID_SUMMARY[7] + ' '.rjust(rjust_width(HEADER_ROW_NAMES_EXPERIMENT_HYBRID_SUMMARY[7]), '~') + str(left_info['correct_parents']) + '\n' + \
               HEADER_ROW_NAMES_EXPERIMENT_HYBRID_SUMMARY[8] + ' '.rjust(rjust_width(HEADER_ROW_NAMES_EXPERIMENT_HYBRID_SUMMARY[8]), '~') + str(percent_of(left_info['correct_parents'], hyb_pep_count)) + '\n' + \
               HEADER_ROW_NAMES_EXPERIMENT_HYBRID_SUMMARY[9] + ' '.rjust(rjust_width(HEADER_ROW_NAMES_EXPERIMENT_HYBRID_SUMMARY[9]), '~') + str(left_info['correct_starting_position']) + '\n' + \
               HEADER_ROW_NAMES_EXPERIMENT_HYBRID_SUMMARY[10] + ' '.rjust(rjust_width(HEADER_ROW_NAMES_EXPERIMENT_HYBRID_SUMMARY[10]), '~') + str(percent_of(left_info['correct_starting_position'], hyb_pep_count)) + '\n' + \
               HEADER_ROW_NAMES_EXPERIMENT_HYBRID_SUMMARY[11] + ' '.rjust(rjust_width(HEADER_ROW_NAMES_EXPERIMENT_HYBRID_SUMMARY[11]), '~') + str(left_info['correct_length']) + '\n' + \
               HEADER_ROW_NAMES_EXPERIMENT_HYBRID_SUMMARY[12] + ' '.rjust(rjust_width(HEADER_ROW_NAMES_EXPERIMENT_HYBRID_SUMMARY[12]), '~') + str(percent_of(left_info['correct_length'], hyb_pep_count)) + '\n' 
    
    summary += '\n' + RIGHT_PARENT_SUMMARY_HEADER + '\n'
    # correct predictions***', (%), near miss****, (%), correct protein, (%), correct starting position, (%), correct length, (%)
    summary += HEADER_ROW_NAMES_EXPERIMENT_HYBRID_SUMMARY[3] + ' '.rjust(rjust_width(HEADER_ROW_NAMES_EXPERIMENT_HYBRID_SUMMARY[3]), '~') + str(right_info['correct_count']) + '\n' + \
               HEADER_ROW_NAMES_EXPERIMENT_HYBRID_SUMMARY[4] + ' '.rjust(rjust_width(HEADER_ROW_NAMES_EXPERIMENT_HYBRID_SUMMARY[4]), '~') + str(percent_of(right_info['correct_count'], hyb_pep_count)) + '\n' + \
               HEADER_ROW_NAMES_EXPERIMENT_HYBRID_SUMMARY[5] + ' '.rjust(rjust_width(HEADER_ROW_NAMES_EXPERIMENT_HYBRID_SUMMARY[5]), '~') + str(right_info['near_miss']) + '\n' + \
               HEADER_ROW_NAMES_EXPERIMENT_HYBRID_SUMMARY[6] + ' '.rjust(rjust_width(HEADER_ROW_NAMES_EXPERIMENT_HYBRID_SUMMARY[6]), '~') + str(percent_of(right_info['near_miss'], hyb_pep_count)) + '\n' + \
               HEADER_ROW_NAMES_EXPERIMENT_HYBRID_SUMMARY[7] + ' '.rjust(rjust_width(HEADER_ROW_NAMES_EXPERIMENT_HYBRID_SUMMARY[7]), '~') + str(right_info['correct_parents']) + '\n' + \
               HEADER_ROW_NAMES_EXPERIMENT_HYBRID_SUMMARY[8] + ' '.rjust(rjust_width(HEADER_ROW_NAMES_EXPERIMENT_HYBRID_SUMMARY[8]), '~') + str(percent_of(right_info['correct_parents'], hyb_pep_count)) + '\n' + \
               HEADER_ROW_NAMES_EXPERIMENT_HYBRID_SUMMARY[9] + ' '.rjust(rjust_width(HEADER_ROW_NAMES_EXPERIMENT_HYBRID_SUMMARY[9]), '~') + str(right_info['correct_starting_position']) + '\n' + \
               HEADER_ROW_NAMES_EXPERIMENT_HYBRID_SUMMARY[10] + ' '.rjust(rjust_width(HEADER_ROW_NAMES_EXPERIMENT_HYBRID_SUMMARY[10]), '~') + str(percent_of(right_info['correct_starting_position'], hyb_pep_count)) + '\n' + \
               HEADER_ROW_NAMES_EXPERIMENT_HYBRID_SUMMARY[11] + ' '.rjust(rjust_width(HEADER_ROW_NAMES_EXPERIMENT_HYBRID_SUMMARY[11]), '~') + str(right_info['correct_length']) + '\n' + \
               HEADER_ROW_NAMES_EXPERIMENT_HYBRID_SUMMARY[12] + ' '.rjust(rjust_width(HEADER_ROW_NAMES_EXPERIMENT_HYBRID_SUMMARY[12]), '~') + str(percent_of(right_info['correct_length'], hyb_pep_count)) + '\n' 
    
    summary_dict['hybrid_peptide_summary'] = {
        'count': hyb_pep_count,
        'totally_correct': total_correct_count, 
        'totally_correct_percent': total_correct_count_p,
        'left_parent': left_info,
        'right_parent': right_info
    }
    
    return summary, summary_dict
    

def __prediction_summary(exp: dict, summary_dict: dict) -> (str, dict):
    '''
    Generate a summary of the predictions done
    
    Inputs:
        exp:     dictionary with the experiment data
        summary_dict: dictionary to save summary information in
    Outputs:
        summary, summary_dict 
            summary:    string of the summary
            summary_dict: same dictionary passed in with new information in it
    '''
    summary = ''
    non_hyb_sum, summary_dict = __non_hybrid_peptide_summary(exp, summary_dict)
    hyb_sum, summary_dict = __hybrid_peptide_summary(exp, summary_dict)
    
    summary += non_hyb_sum + hyb_sum
    
    return summary, summary_dict

######################################################################
#                     END PRIVATE FUNCTIONS
######################################################################

def make_summary(exp: dict, output_dir='./') -> None:
    '''
    Creates a summary txt file for the experiment

    Inputs:
        exp:        dictionary of the experiment results
    kwargs:
        output_dir: string name of directory to save file to. Default='./'
    Outputs:
        None
    '''
    output_dir = __make_valid_dir_string(output_dir)
    __make_dir(output_dir)

    summary_dict = {}
    
    protein_info = exp['experiment_info']['proteins']
    peptide_info = exp['experiment_info']['peptides']

    # start with experiment summary
    # proteins
    # peptides
    header, summary_dict = __summary_header(protein_info, peptide_info, summary_dict)
    
    # crunch the numbers 
    summary, summary_dict = __prediction_summary(exp, summary_dict)
    
    print('Saving summary...')
    with open(output_dir + SUMMARY_FILE_NAME, 'w') as o:
        o.write(header + summary)

    JSON.save_dict(output_dir + SUMMARY_JSON_FILE_NAME, summary_dict)
    print('Done.')