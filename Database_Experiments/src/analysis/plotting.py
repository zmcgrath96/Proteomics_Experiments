from matplotlib import pyplot as plt 
import json
import itertools
import os
import numpy as np
from copy import deepcopy
from sequence_generation.digest import load_digest
from utils import __get_related_files, __make_dir, __make_valid_dir_string, __gzip_dir
from analysis.analysis_utils import get_top_n, __get_argmax_max
from analysis.write_output import write_raw_json, write_summary
from analysis.score_utils import __align_scan_pos, __get_scores_scan_pos_label, __pad_scores

####################################################
#               CONSTANTS
####################################################
plot_colors = [
    'r', 'b', 'g', 'y', 'k', 'c', 'm'
]
plot_markers = [
    ' ', '--', '-.', ':'
]
agg_funcs = ['product', 'sum', 'z_score_sum']
all_line_types = [''.join(x) for x in list(itertools.product(plot_colors, plot_markers))]
all_line_types.sort(key=lambda x: x[1])

cwd = os.path.dirname(os.path.realpath(__file__))
save_fig_prefix = 'figure_{}'

json_header = 'experiment_info'
json_header_prots = 'proteins'
json_header_peps = 'peptides'
json_exp = 'experiment'
json_analysis =  'analysis'
json_rankings = 'ranks'
json_rankings_rankings = 'ranks'
json_seq_len = 'sequence_length'

hybrid_prefix = 'HYBRID'
####################################################
#              END CONSTANTS
####################################################

####################################################
#                 GLOBALS
####################################################
save_fig_count = 0
####################################################
#               END GLOBALS
####################################################

####################################################
#              "PRIVATE" FUNCTIONS
####################################################

'''__plot_subsequence

DESC:
    generate and save a subsequence against the aggregation scores
PARAMS:
    aggs: dictionary of list of aggregation scores. Entry names used for legend
OPTIONAL:
    title: string to label this plot. Also saving name for the plot. Default=''
    save_dir: string name of directory to save all output under. Default =./
    show_graph: bool whether or not to show the graph. Default=False
    peaks: list of dictionaries of the form
            [{
                'protein_name': str,
                'position': int,
                'score': float
            }]
            these are the points of interest to mark
    save_raw_json: bool save the raw json data for each graph. Default=False
    compress: bool to compress the directory. Default=True
    hide_hybrids: bool whether or not to hide hybrid proteins when plotting. Default=True
'''
def __plot_subsequence(aggs, title='', save_dir='./', show_graph=False, agg_func='sum', peaks=None, sequence_info=None, save_raw_json=False, compress=True, hide_hybrids=True):
    save_dir = __make_valid_dir_string(save_dir)
    __make_dir(save_dir)

    plt.figure(figsize=(10, 7))
    i = 0
    for agg in aggs:
        score = aggs[agg]
        if len(score) == 0:
            continue
        if hide_hybrids and hybrid_prefix in agg:
            continue
        plt.plot([j for j in range(len(score))], score, all_line_types[i].strip(), label=str(agg))
        i += 1
    
    if peaks is not None and type(peaks[0]) is dict:
        if hide_hybrids:
            peaks = [p for p in peaks if not hybrid_prefix in p['protein_name']]
        if len(peaks) > 0:
            label_peaks_pos = [x['position'] for x in peaks]
            label_peaks_h = [x['score'] for x in peaks]
            max_peak_value = max([x['score'] for x in peaks])
            max_peak_pos = None
            for x in peaks:
                if x['score'] == max_peak_value:
                    max_peak_pos = x['position']
            plt.plot(label_peaks_pos, label_peaks_h, 'x')
            plt.text(.1, .975, 'max peak position: {}'.format(max_peak_pos), transform=plt.gcf().transFigure)
            plt.text(.1, .96, 'max peak value: {}'.format(max_peak_value),transform=plt.gcf().transFigure)
    if type(sequence_info) is dict:
        plt.text(.5, .975, 'sequence: {}'.format(sequence_info['peptide_sequence']), transform=plt.gcf().transFigure)
        plt.text(.5, .95, 'actual starting position: {}'.format(sequence_info['start_index']), transform=plt.gcf().transFigure)
        plt.text(.5, .925, 'actual parent protein: {}'.format(sequence_info['parent_name']), transform=plt.gcf().transFigure)

    plt.xlabel('subsequence start position')
    plt.ylabel('{} of k-mer scores'.format(agg_func))
    plt.legend()
    plt.title(title)
    plt.savefig(save_dir + title)
    show_graph and plt.show()
    save_raw_json and write_raw_json(save_dir + title, aggs)
    plt.close()
    compress and __gzip_dir(save_dir)

'''plot_subsequence_vs_protein

DESC:
    generate and save a subsequence vs a protein graph
PARAMS:
    k_mers: dictionary of lists of floats to plot. entry names used for legend
OPTIONAL:
    title: string to label this plot. Also saving name for the plot. Default=''
    save_dir: string name of directory to save all output under. Default=./
    show_graph: bool whether or not to show the graph. Default=False
    save_raw_json: bool whether or not to save to raw json data. Default=False
'''
def __plot_subsequence_vs_protein(k_mers, title='', save_dir='./', show_graph=False, save_raw_json=False):
    save_dir = __make_valid_dir_string(save_dir) + title + '/'
    __make_dir(save_dir)

    plt.figure(figsize=(10, 7))
    i  = 0
    for mer in k_mers:
        score = k_mers[mer]
        plt.plot([j for j in range(len(score))], score, all_line_types[i].strip(), label=str(mer))
        i += 1

    plt.xlabel('subsequence start position')
    plt.ylabel('k-mer scores')
    plt.title(title)
    plt.legend()
    plt.savefig(save_dir + title)
    show_graph and plt.show()
    save_raw_json and write_raw_json(save_dir + title, k_mers)
    plt.close()

'''__find_agg_score

DESC:
    find the aggregation
PARAMS:
    scores: dict of list of floats from the analysis
OPTIONAL:
    agg_func: str name of the aggregation function. If none is given, looks for known aggregation functions. Default=''
RETURNS:
    list of floats of the aggregation
'''
def __find_agg_score(scores, agg_func=''):
    known_agg_funcs = [str(x).lower() for x in agg_funcs]

    for key in scores:
        if agg_func == '' or agg_func is None:
            if str(key).lower() in known_agg_funcs:
                return scores[key]
        else:
            if str(agg_func).lower() in str(key).lower():
                return scores[key]

'''__collect_k_rankings

DESC:
    Collect the score rankings for each peptide for each k
PARAMS:
    exp: dictionary experiment summary dictionary
RETURNS:
    dictionary of the form
    {
        k: {
            l_0: [], 
            l_1: [],
            ...
        }
    }
'''
def __collect_k_rankings(exp):
    rankings = {}
    for _, pep in exp[json_exp].items():
        l = pep[json_analysis][json_rankings][json_seq_len]
        for k, rank in pep[json_analysis][json_rankings][json_rankings_rankings].items():
            if not k in rankings:
                rankings[k] = {}
            if not l in rankings[k]:
                rankings[k][l] = []
            rankings[k][l] += [rank]

    return rankings

'''__plot_score_rankings

DESC:
    Generates violin plots of ranking history
PARAMS:
    exp: dictionary experiment summary dictionary
OPTIONAL:
    save_dir: str path to directory to save images too. Default=./
    show_all: bool to show or not all plots. Default=False
RETURNS:
'''
def __plot_score_rankings(exp, save_dir='./', show_all=False):
    # go through each peptide and collect stats
    # x axis: length of the sequnce
    # y axis: violin plot of distribution
    # do it for each k or aggregation
    save_dir = __make_valid_dir_string(save_dir)
    __make_dir(save_dir)

    rs = __collect_k_rankings(exp)
    for k in rs:
        pos = []
        data =[]
        for t in rs[k]:
            pos.append(int(t)) 
            # clean up. Remove None values from lists
            d = [x for x in rs[k][t] if x is not None]
            data.append(d)
        
        if not isinstance(data, list) or not all(isinstance(x, list) for x in data) or not all(isinstance(x, int) for y in data for x in y): 
            print('ERROR: data is not in correct format ( [[int]] ). \ndata: {}'.format(data))
            continue
        if not isinstance(pos, list) or not all(isinstance(x, int) for x in pos): 
            print('ERROR: pos is not in correnct format ( [int] ). \npos: {}'.format(pos))
            continue
   
        plt.violinplot(data, pos)
        plt.title(k)
        plt.xlabel('subsequence length')
        plt.ylabel('score distribution')
        show_all and plt.show()
        plt.savefig(save_dir + str(k))
        plt.close()

#####################################################
#           END "PRIVATE" FUNCTIONS
#####################################################

'''plot_experiment

DESC:
    Plot the score of sequences against the starting position for an experiment
PARAMS:
    exp: dict the json used to save all information
OPTIONAL:
    agg_func: string aggregate function to use. Default=sum
    saving_dir: string path to directory to save figures under. Default=./
    show_all: bool whether or not to show all graphs. Default=False
    use_top_n: bool wether or not to report only the top n scores for a subsequence. Default=False
    n: int the top n scores to keep for a subsequence if use_top_n is True. Default=5
    measure: string the measuring function for determining the top n scores
    compress: bool whether or not to compress output files. If True, all subsequnce plots are compressed. Default=True
    hide_hybrids: bool whether or not to hide hybrid proteins in plots. Default=True
RETURNS: 
    None
'''
def plot_experiment(exp, agg_func='sum', show_all=False, saving_dir='./', use_top_n=False, n=5, measure='average', compress=True, hide_hybrids=True):
    '''
    1. plot the kmer scores
    2. plot the aggregation
    2. plot the ranking history
    '''

    #create the saving directory
    saving_dir = __make_valid_dir_string(saving_dir)
    __make_dir(saving_dir)

    print('\nGenerating plots...')

    # Plot the kmer scores and score aggregations
    # TODO: make this into its own function for the sake of cleanliness
    print('Generating peptide score plots...')
    header_peps = exp[json_header][json_header_peps]

    num_subsequences = len(exp[json_exp])
    subsequence_counter = 0
    for peptide, prots in exp[json_exp].items():
        print('Plotting subsequence {}/{} [{}%]\r'.format(subsequence_counter, num_subsequences, int( (float (subsequence_counter) / float (num_subsequences) ) * 100 )), end="")
        subsequence_counter += 1
        # generate plots for all the proteins against the peptide
        agg_scores = {}
        pep_saving_dir = __make_valid_dir_string(saving_dir + 'subsequence_plots/' + peptide)
        __make_dir(pep_saving_dir)
        for prot, k_mers in prots.items():
            if prot == 'analysis' or prot == 'ranks':
                continue
            plot_title = '{} vs {}'.format(peptide, prot)
            pep_prot_saving_dir = __make_valid_dir_string(pep_saving_dir + prot)
            __make_dir(pep_prot_saving_dir)
            __plot_subsequence_vs_protein(k_mers, title=plot_title, save_dir=pep_prot_saving_dir, show_graph=show_all)
            agg_scores[prot] = __find_agg_score(k_mers, agg_func)
        info = None
        for pep in header_peps:
            if pep['peptide_name'] == peptide:
                info = pep 
                break
        __plot_subsequence(agg_scores, title=str(peptide), save_dir=pep_saving_dir, show_graph=show_all, agg_func=agg_func, peaks=prots['analysis']['predicted_parents'], sequence_info=info, compress=compress, hide_hybrids=hide_hybrids)
    print('Finished.')

    # Plot the ranking of the corect 
    print('Generating score ranking plots...')
    __plot_score_rankings(exp, save_dir=saving_dir + 'ranking_plots/', show_all=show_all)

    print('Finished.')