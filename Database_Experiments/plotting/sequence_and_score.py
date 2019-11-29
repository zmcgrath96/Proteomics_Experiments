from matplotlib import pyplot as plot 
import pandas as pd
import json
import itertools
import os
import numpy as np

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
digests = None
cwd = os.path.dirname(os.path.realpath(__file__))
save_fig_count = 0
save_fig_prefix = 'figure_{}'
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
    try:
        df = pd.read_csv(file, '\t', header=0)
        df = df.sort_values(SCORE_COL, ascending=False)
        df = df.drop_duplicates(subset=SCAN_NO_COL)
        df = df.sort_values(SCAN_NO_COL)
        df = df[df[FILE_NAME_COL].str.contains(search_substring)] if search_substring is not None and search_substring != '' else df

        if not len(df[FILE_NAME_COL]) > 0:
            return [], [], ''
        aligned_scores, _ = __align_scan_pos(list(df[SCORE_COL]), list(df[SCAN_NO_COL]))
        return aligned_scores, [], ''
    except Exception:
        print('could not open file: {}'.format(file))
        return [], [], ''

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

'''__get_related_files

DESC:
    Given some substring, return all files with that substring
PARAMS:
    files: list of strings of file names
    sub: string substring to find in the files
OPTIONAL: 
    not_sub: string a substring that if found in file name don't add
'''
def __get_related_files(files, sub, not_sub=None):
    if not_sub is not None and not_sub != '':
        return [x for x in files if sub in x and not_sub not in x]
    return [x for x in files if sub in x]

'''__load_digest

DESC:
    load the digests from the digestion tsv
PARAMS:
    digest_file: string file path for the digestion tsv
RETURNS:
    None
'''
def __load_digest(digest_file):
    global digests
    if digests is not None:
        return 
    digests = {}
    with open(digest_file, 'r') as o:
        for i, line in enumerate(o):
            if i == 0:
                continue
            l = line.split('\t')
            name = 'peptide_' + str(i-1)
            entry = {'peptide_sequence': l[0], 'parent_name': l[1], 'parent_sequence': l[2], 'start_index': int(l[3])}
            digests[name] = entry 

'''__get_peptide
DESC:
    given a name, get the entry associated with it
PARAMS:
    peptide_name: string name of a peptide entry (peptide_{integer})
RETURNS:
    dictionary entry if it exists else None
'''
def __get_peptide(peptide_name, saving_dir='./'):
    global digests
    saving_dir = saving_dir + '/' if not saving_dir[-1] == '/' else saving_dir
    if not digests:
        file_name =  saving_dir + 'digestion.tsv'
        __load_digest(file_name)
    return digests[peptide_name] if peptide_name in digests else None

'''__parse_output_name

DESC:
    get a usable label for a comparision 
PARAMS: 
    crux_output: a string that will use the parent dir of the output as the label
RETURNS:
    a string with the title 
'''
def __parse_output_name(crux_output):
    return crux_output.split('/')[-2]

####################################################
#              END UTILITY FUNCTIONS
####################################################

'''plot_subsequence

DESC:
    generate and save all graphs for a subsequence
PARAMS:
    files: list of strings each entry is an output crux file who's contents are related
OPTIONAL:
    title: string to label this plot. Default=''
    save_dir: string name of directory to save all output under. Default=./
    aggregate: string which type of aggregation of scores to perfrom. Default=sum
    show_graph: bool whether or not to show the graph
'''
def plot_subsequence(files, title='', save_dir='./', aggregate='sum', show_graph=False):
    global save_fig_count, save_fig_prefix
    if title is None or title == '':
        title = save_fig_prefix.format(save_fig_count)
        save_fig_count += 1

    # create the saving directory if it doesn't exits
    if not os.path.exists(save_dir): 
        os.makedirs(save_dir)

    plot.figure(figsize=(10,7))
    total_score, all_scores = score_vs_position(files, aggregate=aggregate)
    i = 0
    for lab, scores in all_scores.items():
        this_label = str([int(i) for i in str(lab).replace('_', ' ').split() if i.isdigit()][0]) + '-mer'
        plot.plot([j for j in range(len(scores))], scores, all_line_types[i].strip(), label=this_label)
        i += 1
    plot.plot([j for j in range(len(total_score))], total_score, all_line_types[len(all_scores)], label='{} of scores'.format(aggregate))
    plot.legend()
    plot.title(title)
    
    plot.xlabel('k-mer starting position')
    plot.ylabel('k-mer score')
    plot.savefig(save_dir + title)
    show_graph and plot.show()
    return total_score

'''score_vs_position

DESC:
    Plot the peptide score vs position in the sequnce
PARAMS:
    files -- list of strings to crux output files
RETRUNS:
    total_scores: aggregate of all scores into one list
    all_scores: a dictionary of lists where the entry key is the label for the list
'''
def score_vs_position(files, aggregate='sum', search_substring=''):
    total_score = []
    all_scores = {}
    for score_file in files:
        # get scores, scan numbers, and labels
        scores, _, _ = __get_scores_scan_pos_label(score_file, search_substring)
        # if there isn't any data, continue
        if not len(scores) > 0:
            continue
        this_label = __parse_output_name(score_file)
        all_scores[this_label] = scores
        total_score, scores = __pad_scores(total_score, scores)
        for j in range(len(scores)):
            if aggregate == 'product' or 'product' in aggregate.lower():
                if scores[j] <= 0:
                    continue
                total_score[j] *= scores[j]
            else:
                total_score[j] += scores[j]

    return total_score, all_scores

'''plot_experiment

DESC:
    Plot the score of sequences against the starting position for an experiment
PARAMS:
    experiment: a string to determine which experiemnt to plot. Defaults flipped
    files: a list of strings of file paths for results
    sequences_json: string for path to json that holds parent info
'''

def plot_experiment(experiment, files, protein_names, subsequence_prefix, num_subsequences, hybrid_prefix, agg_func='sum', show_all=False, saving_dir='./'):
    saving_dir = saving_dir + '/' if saving_dir[-1] != '/' else saving_dir

    #create the saving directory
    if not os.path.exists(saving_dir): 
        os.makedirs(saving_dir)

    if 'fractionated' in str(experiment).lower():
        pass

    else:
        print('\nGenerating plots...')
        # hybrid first
        hybrid_related = __get_related_files(files, hybrid_prefix)
        total_scores = {}
        for protein_name in protein_names:
            prot_with_hybrid = __get_related_files(hybrid_related, protein_name)
            total_scores[protein_name] = plot_subsequence(prot_with_hybrid, title='hybrid vs {}'.format(protein_name), aggregate=agg_func, save_dir=saving_dir + 'hybrid_{}/'.format(protein_name), show_graph=show_all)
        for i, title in enumerate(total_scores):
            plot.plot([j for j in range(len(total_scores[title]))], total_scores[title], all_line_types[i].strip(), label=title)
        
        plot.legend()
        plot.title('hybrid sequence')
        plot.xlabel('k-mer starting position')
        plot.ylabel('{} of k-mer scores'.format(agg_func))
        plot.figure(figsize=(10,7))
        show_all and plot.show()
        plot.savefig(saving_dir + 'hybrid')
        #clear up variables
        hybrid_related = None 
        total_scores = None 

        # do the rest
        peptide_names = ['{}_{}'.format(subsequence_prefix, x) for x in range(num_subsequences)]
        for i, peptide_name in enumerate(peptide_names):
            print('Generating graphs for peptide sequences {}/{}[{}%]'.format(i+1, len(peptide_names), int(((i+1)/len(peptide_names)*100))))
            peptide_entry = __get_peptide(peptide_name, saving_dir=saving_dir)
            pep_related = __get_related_files(files, peptide_name)
            total_scores = {}
            max_score_position = 0
            max_score = 0
            for protein_name in protein_names:
                prot_with_pep = __get_related_files(pep_related, protein_name)
                total_scores[protein_name] = plot_subsequence(prot_with_pep, title='{} vs {}'.format(peptide_name, protein_name), aggregate=agg_func, save_dir=saving_dir+ '{}_{}/'.format(peptide_name, protein_name), show_graph=show_all)
            for i, title in enumerate(total_scores):
                if max(total_scores[title]) > max_score:
                    max_score = max(total_scores[title])
                    max_score_position = np.argmax(total_scores[title])

                plot.plot([j for j in range(len(total_scores[title]))], total_scores[title], all_line_types[i].strip(), label=title)
            
            plot.legend()
            plot.title(peptide_name)
            plot.figtext(0.55, 0.98, 'peptide sequence: {}'.format(peptide_entry['peptide_sequence']))
            plot.figtext(0.55, 0.96, 'parent protein: {}'.format(peptide_entry['parent_name']))
            plot.figtext(0.55, 0.94, 'peptide starting position in parent: {}'.format(peptide_entry['start_index']))
            plot.figtext(0.55, 0.92, 'maximum scoring position: {}'.format(max_score_position))
            plot.xlabel('k-mer starting position')
            plot.ylabel('{} of k-mer scores'.format(agg_func))
            plot.figure(figsize=(10,7))
            show_all and plot.show()
            plot.savefig(saving_dir + peptide_name)

            #clear up variables 
            total_scores = None
            pep_related = None 
        print('Finished')

        


