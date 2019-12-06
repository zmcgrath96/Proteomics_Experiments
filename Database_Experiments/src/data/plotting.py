from matplotlib import pyplot as plt 
import json
import itertools
import os
import numpy as np
from copy import deepcopy
from sequences.digest import load_digest
from utils import __get_related_files, __make_dir, __make_valid_dir_string
from data.analysis import get_top_n, __get_argmax_max
from data.write_output import write_raw_json, write_summary
from data.score_utils import __align_scan_pos, __get_scores_scan_pos_label, __pad_scores
from data.aggregations import __sum, __product, __z_score_sum

####################################################
#               CONSTANTS
####################################################
plot_colors = [
    'r', 'b', 'g', 'y', 'k', 'c', 'm'
]
plot_markers = [
    ' ', '--', '-.', ':'
]
all_line_types = [''.join(x) for x in list(itertools.product(plot_colors, plot_markers))]
all_line_types.sort(key=lambda x: x[1])

cwd = os.path.dirname(os.path.realpath(__file__))
save_fig_prefix = 'figure_{}'

json_header = 'header'
json_header_prots = 'proteins'
json_header_peps = 'peptides'
json_exp = 'experiment'
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
'''
def __plot_subsequence(aggs, title='', save_dir='./', show_graph=False, agg_func='sum'):
    save_dir = __make_valid_dir_string(save_dir)
    __make_dir(save_dir)

    plt.figure(figsize=(10, 7))
    i = 0
    for agg in aggs:
        score = aggs[agg]
        if len(score) == 0:
            continue
        plt.plot([j for j in range(len(score))], score, all_line_types[i].strip(), label=str(agg))
        i += 1
    
    plt.xlabel('subsequence start position')
    plt.ylabel('{} of k-mer scores'.format(agg_func))
    plt.legend()
    plt.title(title)
    show_graph and plt.show()
    plt.savefig(save_dir + title)
    plt.close()
    write_raw_json(save_dir + title, aggs)

'''plot_subsequence_vs_protein

DESC:
    generate and save a subsequence vs a protein graph
PARAMS:
    k_mers: dictionary of lists of floats to plot. entry names used for legend
OPTIONAL:
    title: string to label this plot. Also saving name for the plot. Default=''
    save_dir: string name of directory to save all output under. Default=./
    aggregate: string which type of aggregation of scores to perfrom. Default=sum
    show_graph: bool whether or not to show the graph. Default=False
'''
def __plot_subsequence_vs_protein(k_mers, title='', save_dir='./', aggregate='sum', show_graph=False):
    agg = __z_score_sum if 'z_score_sum' in aggregate.lower() else (__product if 'product' in aggregate.lower() else __sum)
    save_dir = __make_valid_dir_string(save_dir) + title + '/'
    __make_dir(save_dir)

    plt.figure(figsize=(10, 7))
    i  = 0
    all_scores = []
    for mer in k_mers:
        score = k_mers[mer]
        plt.plot([j for j in range(len(score))], score, all_line_types[i].strip(), label=str(mer))
        i += 1
        all_scores.append(score)
    agg_scores = agg(all_scores)
    plt.plot([j for j in range(len(agg_scores))], agg_scores, all_line_types[i], label='{} of scores'.format(aggregate))
    plt.xlabel('subsequence start position')
    plt.ylabel('k-mer scores')
    plt.title(title)
    plt.legend()
    show_graph and plt.show()
    plt.savefig(save_dir + title)
    plt.close()
    write_raw_json(save_dir + title, k_mers)
    return agg_scores

#####################################################
#           END "PRIVATE" FUNCTIONS
#####################################################

'''plot_experiment

DESC:
    Plot the score of sequences against the starting position for an experiment
PARAMS:
    experiment: string name of experiment to run
    experiment_json_file: string path to the experiment json file saved to disk
OPTIONAL:
    agg_func: string aggregate function to use. Default=sum
    saving_dir: string path to directory to save figures under. Default=./
    show_all: bool whether or not to show all graphs. Default=False
    use_top_n: bool wether or not to report only the top n scores for a subsequence. Default=False
    n: int the top n scores to keep for a subsequence if use_top_n is True. Default=5
    measure: string the measuring function for determining the top n scores
RETURNS: 
    None
'''
def plot_experiment(experiment, experiment_json_file, agg_func='sum', show_all=False, saving_dir='./', use_top_n=False, n=5, measure='average'):
    #create the saving directory
    saving_dir = __make_valid_dir_string(saving_dir)
    __make_dir(saving_dir)
    exp = {}
    with open(experiment_json_file, 'r') as exp_file:
        exp = json.load(exp_file)

    if 'fractionated' in str(experiment).lower():
        pass

    else:
        print('\nGenerating plots...')
        for peptide, prots in exp[json_exp].items():
            # generate plots for all the proteins against the peptide
            agg_scores = {}
            pep_saving_dir = __make_valid_dir_string(saving_dir + peptide)
            agg_scores[peptide] = {}
            __make_dir(pep_saving_dir)
            for prot, k_mers in prots.items():
                plot_title = '{} vs {}'.format(peptide, prot)
                pep_prot_saving_dir = __make_valid_dir_string(pep_saving_dir + prot)
                __make_dir(pep_prot_saving_dir)
                agg_scores[prot] = __plot_subsequence_vs_protein(k_mers, title=plot_title, save_dir=pep_prot_saving_dir, aggregate=agg_func, show_graph=show_all)
            __plot_subsequence(agg_scores, title=str(peptide), save_dir=pep_saving_dir, show_graph=show_all, agg_func=agg_func)
        print('Finished')