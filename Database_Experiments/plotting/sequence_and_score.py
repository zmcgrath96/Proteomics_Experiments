from matplotlib import pyplot
import pandas as pd

FILE_NAME_COL = 'file'
SCORE_COL = 'xcorr score'
SCAN_NO_COL = 'scan'

plot_colors = [
    'r', 'b', 'g', 'y', 'k', 'c', 'm'
]

'''score_vs_position

DESC:
    Plot the peptide score vs position in the sequnce
PARAMS:
    files -- list of strings to crux output files
RETRUNS:
    None
'''
def score_vs_position(files):
    if len(files) > len(plot_colors):
        for score_file in files:
            df = pd.read_csv(score_file, '\t', header=0)
            df = df.sort_values(SCORE_COL, ascending=False)
            df = df.drop_duplicates(subset=SCAN_NO_COL)
            df = df.sort_values(SCAN_NO_COL)
            title = df[FILE_NAME_COL][0].split('/')[-1]
            pyplot.plot(df[SCAN_NO_COL], df[SCORE_COL], 'ro')
            pyplot.xlabel('subsequence start position')
            pyplot.ylabel(SCORE_COL)
            pyplot.title(title)
            pyplot.show()
    else:
        maxes = []
        scan_lens = []
        for i, score_file in enumerate(files):
            df = pd.read_csv(score_file, '\t', header=0)
            df = df.sort_values(SCORE_COL, ascending=False)
            df = df.drop_duplicates(subset=SCAN_NO_COL)
            df = df.sort_values(SCAN_NO_COL)
            if not len(df[FILE_NAME_COL]) > 0:
                continue
            maxes.append(max(df[SCORE_COL]))
            scan_lens.append(max(df[SCAN_NO_COL]))
            this_label = df[FILE_NAME_COL][0].split('/')[-1] 
            pyplot.plot(df[SCAN_NO_COL], df[SCORE_COL], plot_colors[i], label=this_label)
        pyplot.legend()
        pyplot.show()

