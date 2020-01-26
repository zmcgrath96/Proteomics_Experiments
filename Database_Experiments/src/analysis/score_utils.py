import pandas as pd

###############################################
#               CONSTANTS 
###############################################
FILE_NAME_COL = 'file'
SCORE_COL = 'score'
SCAN_NO_COL = 'scan_no'
###############################################
#             END CONSTANTS
###############################################

'''__align_scan_pos

DESC:
    Align scores list to the correspondeing position from the scan numbers
PARAMS:
    scores: list of floats of scores
    scan_nos: list of ints of positions where the scores should go
RETURNS:
    list, list: aligned scores, aligned scan positions
'''
def __align_scan_pos(scores, scan_nos):
    aligned_scores = [0 for _ in range(max(scan_nos) + 1)]
    aligned_scan_nos = [i for i in range(max(scan_nos) + 1)]
    for i, insert_index in enumerate(scan_nos):
        # NOTE: for the mzML files, the first scan starts at 1, not at 0 so minus 1
        # NOTE 2: when i updated to write own scoring the number was off again so removing th e-1
        aligned_scores[insert_index] = scores[i]

    return aligned_scores, aligned_scan_nos

'''__get_scores_scan_pos_label

DESC:
    Extract the scores and the scan positions from the file
PARAMS:
    file: a string for the filepath to extract results from
OPTIONAL:
    search_substring: a string to search through filenames to limit the search. Defaults to empty
RETUNS:
    list, list: lists of the scores, scan numbers
'''
def __get_scores_scan_pos_label(file, search_substring=''):
    try:
        sep = '\t' if '.tsv' in file else ','
        df = pd.read_csv(file, sep, header=0)
        df = df.sort_values(SCORE_COL, ascending=False)
        df = df.drop_duplicates(subset=SCAN_NO_COL)
        df = df.sort_values(SCAN_NO_COL)
        #df = df[df[FILE_NAME_COL].str.contains(search_substring)] if search_substring is not None and search_substring != '' else df

        # if not len(df[FILE_NAME_COL]) > 0:
        #     return [], [], ''
        aligned_scores, _ = __align_scan_pos(list(df[SCORE_COL]), list(df[SCAN_NO_COL]))
        return aligned_scores, [], ''
    except Exception:
        print('could not open file: {}'.format(file))
        return [], [], ''

'''__pad_scores

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
def __pad_scores(score_l1, score_l2, padding=0):
    score_l1 = list(score_l1)
    score_l2 = list(score_l2)
    diff = len(score_l1) - len(score_l2)
    if diff == 0:
        return score_l1, score_l2
    elif diff > 0:
        return score_l1, score_l2 + [padding for _ in range(abs(diff))]
    else:
        return score_l1 + [padding for _ in range(abs(diff))], score_l2