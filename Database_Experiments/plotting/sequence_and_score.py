from matplotlib import pyplot
import pandas as pd
import json
import itertools

out_fig_count = 0
####################################################
#               CONSTANTS
####################################################
FILE_NAME_COL = 'file'
SCORE_COL = 'xcorr score'
SCAN_NO_COL = 'scan'

plot_colors = [
    'r', 'b', 'g', 'y', 'k', 'c', 'm'
]
plot_markers = [
    ' ', '--', '-.', ':'
]
all_line_types = [''.join(x) for x in list(itertools.product(plot_colors, plot_markers))]
all_line_types.sort(key=lambda x: x[1])
####################################################
#              END CONSTANTS
####################################################

####################################################
#               UTILITY FUNCTIONS
####################################################

'''__get_scores_and_scan_pos

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
    df = pd.read_csv(file, '\t', header=0)
    df = df.sort_values(SCORE_COL, ascending=False)
    df = df.drop_duplicates(subset=SCAN_NO_COL)
    df = df.sort_values(SCAN_NO_COL)
    df = df[df[FILE_NAME_COL].str.contains(search_substring)] if search_substring is not None and search_substring != '' else df

    if not len(df[FILE_NAME_COL]) > 0:
        return [], [], ''
    aligned_scores, aligned_scan_pos = __align_scan_pos(list(df[SCORE_COL]), list(df[SCAN_NO_COL]))
    return aligned_scores, aligned_scan_pos, str(str(df[FILE_NAME_COL][0].split('/')[-1]).split('_')[-1]).split('.')[0] + '-mer'

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
        aligned_scores[insert_index] = scores[i]

    return aligned_scores, aligned_scan_nos

'''__pad_scores

DESC:
    Makes two lists the same length (filled with 0s)
PARAMS:
    score_l1: list of floats for the first set of scores
    score_l2: list of floats for the second set of scores
RETURNS:
    list, list of modified score_l1, modified score_l2
'''
def __pad_scores(score_l1, score_l2):
    diff = len(score_l1) - len(score_l2)
    if diff == 0:
        return score_l1, score_l2
    elif diff > 0:
        return score_l1, score_l2 + [0 for _ in range(diff)]
    else:
        return score_l1 + [1 for _ in range(abs(diff))], score_l2

####################################################
#              END UTILITY FUNCTIONS
####################################################

'''score_vs_position

DESC:
    Plot the peptide score vs position in the sequnce
PARAMS:
    files -- list of strings to crux output files
RETRUNS:
    None
'''
def score_vs_position(files, search_substring='', title=''):
    i = 0
    total_score = []
    for score_file in files:
        # get scores, scan numbers, and labels
        scores, scan_nos, this_label = __get_scores_scan_pos_label(score_file, search_substring)
        # if there isn't any data, continue
        if not len(scores) > 0:
            continue
        # plot any score we have
        pyplot.plot(scan_nos, scores, all_line_types[i].strip(), label=this_label)
        i += 1
        total_score, scores = __pad_scores(total_score, scores)
        for j in range(len(scores)):
            if scores[j] <= 0:
                continue
            total_score[j] *= scores[j]
    # plot the total sum of the scores
    pyplot.plot([j for j in range(len(total_score))], total_score, all_line_types[i].strip(), label='Product of scores')
    
    global out_fig_count
    title = 'output figure ' + str(out_fig_count) if title is None or title == '' else title
    save_name =  './' + title
    out_fig_count += 1

    pyplot.title(title)
    pyplot.legend()
    pyplot.savefig(save_name)
    pyplot.show()

'''plot_experiment

DESC:
    Plot the score of sequences against the starting position for an experiment
PARAMS:
    experiment: a string to determine which experiemnt to plot. Defaults flipped
    files: a list of strings of file paths for results
    sequences_json: string for path to json that holds parent info
'''

def plot_experiment(experiment, files, sequences_json):
    sequences = None
    with open(sequences_json, 'r') as seqfile:
        sequences = json.load(seqfile)

    if 'fractionated' in str(experiment).lower():
        score_vs_position(files)

    else:
        score_vs_position(files, sequences['parents']['left_parent']['name'], sequences['parents']['left_parent']['name'])
        score_vs_position(files, sequences['parents']['right_parent']['name'], sequences['parents']['right_parent']['name'])
