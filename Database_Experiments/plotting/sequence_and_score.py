from matplotlib import pyplot
import pandas as pd
import json
import itertools

FILE_NAME_COL = 'file'
SCORE_COL = 'xcorr score'
SCAN_NO_COL = 'scan'

plot_colors = [
    'r', 'b', 'g', 'y', 'k', 'c', 'm'
]
plot_markers = [
    ' ', '^', 'o', '+', 'x', '*'
]
all_line_types = [''.join(x) for x in list(itertools.product(plot_colors, plot_markers))]
all_line_types.sort(key=lambda x: x[1])
'''score_vs_position

DESC:
    Plot the peptide score vs position in the sequnce
PARAMS:
    files -- list of strings to crux output files
RETRUNS:
    None
'''
def score_vs_position(files, search_substring=''):
    if len(files) > len(all_line_types):
        for score_file in files:
            df = pd.read_csv(score_file, '\t', header=0)
            df = df.sort_values(SCORE_COL, ascending=False)
            df = df.drop_duplicates(subset=SCAN_NO_COL)
            df = df.sort_values(SCAN_NO_COL)
            df = df[df[FILE_NAME_COL].str.contains(search_substring)] if search_substring is not None and search_substring != '' else df

            title = df[FILE_NAME_COL][0].split('/')[-1]
            pyplot.plot(df[SCAN_NO_COL], df[SCORE_COL], 'ro')
            pyplot.xlabel('subsequence start position')
            pyplot.ylabel(SCORE_COL)
            pyplot.title(title)
            pyplot.show()
    else: 
        i = 0
        for score_file in files:
            df = pd.read_csv(score_file, '\t', header=0)
            df = df.sort_values(SCORE_COL, ascending=False)
            df = df.drop_duplicates(subset=SCAN_NO_COL)
            df = df.sort_values(SCAN_NO_COL)
            df = df[df[FILE_NAME_COL].str.contains(search_substring)] if search_substring is not None and search_substring != '' else df

            if not len(df[FILE_NAME_COL]) > 0:
                continue
            this_label = df[FILE_NAME_COL][0].split('/')[-1] 
            pyplot.plot(df[SCAN_NO_COL], df[SCORE_COL], all_line_types[i].strip(), label=this_label)
            i += 1
        pyplot.legend()
        pyplot.show()

def plot_experiment(experiment, files, sequences_json):
    sequences = None
    with open(sequences_json, 'r') as seqfile:
        sequences = json.load(seqfile)

    if 'fractionated' in str(experiment).lower():
        pass 

    else:
        score_vs_position(files, sequences['parents']['left_parent']['name'])
        score_vs_position(files, sequences['parents']['right_parent']['name'])
