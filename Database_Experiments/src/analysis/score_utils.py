import pandas as pd
import utils

###############################################
#               CONSTANTS 
###############################################
FILE_NAME_COL = 'file'
SCORE_COL = 'score'
SCAN_NO_COL = 'scan_no'

score_funcs = ['custom', 'crux']
column_names = {
    'custom': {
        'file_name': 'file',
        'score': 'score',
        'scan_number': 'scan_no'
    },
    'crux': {
        'file_name': 'file',
        'score': 'xcorr score',
        'scan_number': 'scan'
    }
}
###############################################
#             END CONSTANTS
###############################################

###############################################
#           PRIVATE FUNCTIONS
###############################################
def __get_correct_col_names(col_names):
    d_to_l = lambda d: [x for _, x in d.items()]
    for s_func, d in column_names.items():
        ks = d_to_l(d)
        if set.issubset(set(ks), set(col_names)):
            return s_func
    return 'custom'

###############################################
#           END PRIVATE FUNCTIONS
###############################################


'''align_scan_pos

DESC:
    Align scores list to the correspondeing position from the scan numbers
PARAMS:
    scores: list of floats of scores
    scan_nos: list of ints of positions where the scores should go
RETURNS:
    list, list: aligned scores, aligned scan positions
'''
def align_scan_pos(scores, scan_nos):
    aligned_scores = [0 for _ in range(max(scan_nos) + 1)]
    aligned_scan_nos = [i for i in range(max(scan_nos) + 1)]
    for i, insert_index in enumerate(scan_nos):
        # NOTE: for the mzML files, the first scan starts at 1, not at 0 so minus 1
        # NOTE 2: when i updated to write own scoring the number was off again so removing th e-1
        aligned_scores[insert_index] = scores[i]

    return aligned_scores, aligned_scan_nos

'''get_scores_scan_pos_label

DESC:
    Extract the scores and the scan positions from the file
PARAMS:
    file: a string for the filepath to extract results from
OPTIONAL:
    search_substring: a string to search through filenames to limit the search. Defaults to empty
RETUNS:
    list, list: lists of the scores, scan numbers
'''
def get_scores_scan_pos_label(file, search_substring=''):
    sep = '\t' if '.tsv' in file else ','
    if utils.__is_gzipped(file):
        file = utils.__gunzip(file)

    df = pd.read_csv(file, sep, header=0)
        
    s_func = __get_correct_col_names(df.columns)
    col_names = column_names[s_func]

    df = df.sort_values(col_names['score'], ascending=False)
    df = df.drop_duplicates(subset=col_names['scan_number'])
    df = df.sort_values(col_names['scan_number'])
    
    aligned_scores, _ = align_scan_pos(list(df[col_names['score']]), list(df[col_names['scan_number']]))
    return aligned_scores, [], ''

'''pad_scores

DESC:
    Makes two lists the same length (filled with 0s)
PARAMS:
    score_l1: list of floats for the first set of scores
    score_l2: list of floats for the second set of scores
OPTIONAL:    
    padding: float number to pad the list with. Default=0
RETURNS:
    list, list of modified score_l1, modified score_l2
'''
def pad_scores(score_l1, score_l2, padding=0):
    score_l1 = list(score_l1)
    score_l2 = list(score_l2)
    diff = len(score_l1) - len(score_l2)
    if diff == 0:
        return score_l1, score_l2
    elif diff > 0:
        return score_l1, score_l2 + [padding for _ in range(abs(diff))]
    else:
        return score_l1 + [padding for _ in range(abs(diff))], score_l2