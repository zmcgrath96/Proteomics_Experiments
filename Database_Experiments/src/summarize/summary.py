from numbers import Number
from utils import __make_valid_dir_string, __make_dir
from file_io import JSON

######################################################################
#                   CONSTANTS
######################################################################

DIVIDER_WIDTH = 70
SUMMARY_FILE_NAME = 'summary.txt'
SUMMARY_JSON_FILE_NAME = 'summary.json'
SECTION_DIVIDER = '\n\n' + '=' * DIVIDER_WIDTH + '\n\n\n\n\n'
HEADER_UNDERLINE = '\n' + '-' * DIVIDER_WIDTH + '\n'
PROTEIN_SUMMARY_HEADER = 'EXPERIMENT PROTEIN INFORMATION'
PEPTIDE_SUMMARY_HEADER = 'EXPERIMENT PEPTIDE INFORMATION'
NON_HYBRID_SUMMARY_HEADER = 'NON HYBRID PEPTIDE SUMMARY'
HYBRID_SUMMARY_HEADER = 'HYBRID PEPTIDE SUMMARY'
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
    'total peptides',               # int
    'correct predictions',          # int
    '(%)',                          # float
    'correct parent protein',       # int
    '(%)',                          # float
    'near miss predictions*',       # float
    '(%)',                          # float
    'correct starting position*',   # int
    '(%)',                          # float
    'correct peptide length*',      # int
    '(%)'                           # float
]

__cols = lambda x: ''.join(['{}\t' for _ in range(x)])

HEADER_ROW_PROTEIN_SUMMARY = __cols(len(HEADER_ROW_NAMES_PROTEIN_SUMMARY)).format(*HEADER_ROW_NAMES_PROTEIN_SUMMARY)
HEADER_ROW_PEPTIDE_SUMMARY = __cols(len(HEADER_ROW_NAMES_PEPTIDE_SUMMARY)).format(*HEADER_ROW_NAMES_PEPTIDE_SUMMARY)
HEADER_ROW_EXPERIMENT_NON_HYBRID_SUMMARY = __cols(len(HEADER_ROW_NAMES_EXPERIMENT_NON_HYBRID_SUMMARY)).format(*HEADER_ROW_NAMES_EXPERIMENT_NON_HYBRID_SUMMARY)

NEAR_MISS_CONTEXT = '\n\n* near miss \n \
A near miss means that the predicted peptide sequence identified the correct \n \
parent protein and was ONLY one of the following:\n \
    (1)  Predicted length was off by 1\n \
    (2)  The starting position was off by 1 position\n'

######################################################################
#                         END CONSTANTS
######################################################################

######################################################################
#                       PRIVATE FUNCTIONS
######################################################################

def __pad_and_center(num: Number, s: str) -> str:
    '''
    pad and center a number in a string
    
    Inputs:
        num:        the number to center
        s:          the string to center the number int
    Outputs:
        centered:   string representing the number
    '''
    if len(str(num)) >= len(s):
        return str(num) 
    
    center = ''
    side_len = (len(s) - len(str(num))) // 2 
    center += ' ' * side_len + str(num) + ' ' * side_len
    return center

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

def __summary_header(prot_info: dict, pep_info: dict, summary_dict: dict) -> (str, dict):
    '''
    Create a string to return that has the basic experiment information in it
    
    Inputs:
        prot_info:  dictionary of protein information
        pep_info:   dictionary of peptide information
    Outputs:
        header:     string of the header information
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

def __prediction_summary(exp: dict, summary_dict: dict) -> (str, dict):
    '''
    Generate a summary of the predictions done
    
    Inputs:
        exp:        dictionary with the experiment data
    Outputs:
        summary:    string of the summary
    '''
    #total peptides, correct predictions, (%), near miss predictions (%), correct parent protein, (%), correct starting position, (%), correct peptide length, (%)
    summary = ''
    summary += NON_HYBRID_SUMMARY_HEADER + HEADER_UNDERLINE
    
    # do the number crunching
    peptide_count = len([x for x in exp['experiment_info']['peptides'] if 'hybrid' not in x['peptide_name'].lower()])
    percent_of = lambda x, c: round((float(x)/float(c)) * 100, 2)
    
    corr_c, near_miss, corr_par, corr_start_pos, corr_len = __prediction_stats_non_hyb(exp)
    corr_c_p = percent_of(corr_c, peptide_count)
    near_miss_p = percent_of(near_miss, peptide_count)
    corr_par_p = percent_of(corr_par, peptide_count)
    corr_start_pos_p = percent_of(corr_start_pos, peptide_count)
    corr_len_p = percent_of(corr_len, peptide_count)
    
    rjust_width = lambda s: DIVIDER_WIDTH - len(s)
    
    summary += HEADER_ROW_NAMES_EXPERIMENT_NON_HYBRID_SUMMARY[0] + ' '.rjust(rjust_width(HEADER_ROW_NAMES_EXPERIMENT_NON_HYBRID_SUMMARY[0]), '~') + str(peptide_count) + '\n' +\
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
        'count': peptide_count,
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

######################################################################
#                   END PRIVATE FUNCTIONS
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
    
    with open(output_dir + SUMMARY_FILE_NAME, 'w') as o:
        o.write(header + summary)

    JSON.save_dict(output_dir + SUMMARY_JSON_FILE_NAME, summary_dict)