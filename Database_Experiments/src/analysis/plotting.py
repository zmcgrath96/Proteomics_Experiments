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
from analysis import peptide_plotting, protein_plotting

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

agg_funcs = ['product', 'sum', 'z_score_sum']

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
'''__find_agg_func
'''
def __find_agg_func(scores): 
    known_agg_funcs = [str(x).lower() for x in agg_funcs]
    for key in scores:
        if str(key).lower() in known_agg_funcs:
            return key

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
    agg_func = agg_func if agg_func != '' else __find_agg_func(scores)

    for key in scores:
        if str(agg_func).lower() == str(key).lower():
            return scores[key]

'''__find_pep_starting_pos

'''
def __find_pep_starting_pos(peps, pep_name):
    for pep in peps:
        if pep['peptide_name'] == pep_name:
            return pep['start_index']
    return -1

'''__aggregate_ranks

DESC:
    sort out the ranks for each peptide from the experiment json
PARAMS:
    exp: dictionary that holds the experiment information
RETURNS:
    list of dictionaries of peptides with their ranks
'''
def __aggregate_ranks(exp):
   
    rs = []
    peps_header = exp[json_header][json_header_peps]
    for pep_name, pep in exp[json_exp].items(): 
        pep_rank = pep[json_analysis][json_rankings]
        is_hybrid = True if hybrid_prefix in pep_name else False
        r = {
            'starting_position': __find_pep_starting_pos(peps_header, pep_name),
            'ranking': pep_rank[json_rankings_rankings][__find_agg_func(pep_rank[json_rankings_rankings])], 
            'parent_protein': pep_rank['correct_protein'],
            'sequence_length': pep_rank['sequence_length'],
            'hybrid': is_hybrid
        }
        rs.append(r)
    return rs


#####################################################
#           END "PRIVATE" FUNCTIONS
#####################################################

'''plot_protein_summary

DESC:
    call the protein plotting functions 
PARAMS:
    exp: dict the json used to save all information 
OPTIONAL:
    saving_dir: string path to directory to save figures under. Default=./
    show_all: bool whether or not to show all graphs. Default=False
    compress: bool whether or not to compress output files. If True, all subsequnce plots are compressed. Default=True
RETURNS:
    None
'''
def plot_protein_summary(exp, saving_dir='./', show_all=False, compress=True):
    # clean up proteins (normalize because of the hybrids)
    dirty_prots = exp[json_header][json_header_prots]
    clean_prots = []
    for prot in dirty_prots:
        if 'protein' in prot:
            tp = {
                'name': prot['name'],
                'sequence': prot['protein']
            }
            clean_prots.append(tp)
        else:
            clean_prots.append(prot)
    protein_plotting.prots_pep_pos_rankings(clean_prots, __aggregate_ranks(exp), save_dir=saving_dir, show_all=show_all, compress=compress)
    

'''plot_peptide_scores

DESC:
    Plot the score of sequences against the starting position for an experiment
PARAMS:
    exp: dict the json used to save all information
OPTIONAL:
    agg_func: string aggregate function to use. Default=sum
    saving_dir: string path to directory to save figures under. Default=./
    show_all: bool whether or not to show all graphs. Default=False
    compress: bool whether or not to compress output files. If True, all subsequnce plots are compressed. Default=True
    hide_hybrids: bool whether or not to hide hybrid proteins in plots. Default=True
RETURNS: 
    None
'''
def plot_peptide_scores(exp, agg_func='sum', saving_dir='./', show_all=False, compress=True, hide_hybrids=True):
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
            peptide_plotting.plot_subsequence_vs_protein(k_mers, title=plot_title, save_dir=pep_prot_saving_dir, show_graph=show_all)
            agg_scores[prot] = __find_agg_score(k_mers, agg_func)
        info = None
        for pep in header_peps:
            if pep['peptide_name'] == peptide:
                info = pep 
                break
        peptide_plotting.plot_subsequence(agg_scores, title=str(peptide), save_dir=pep_saving_dir, show_graph=show_all, agg_func=agg_func, peaks=prots['analysis']['predicted_parents'], sequence_info=info, compress=compress, hide_hybrids=hide_hybrids)
    

'''plot_experiment

DESC:
    Generate various plots for the experiment
PARAMS:
    exp: dict the json used to save all information
OPTIONAL:
    agg_func: string aggregate function to use. Default=sum
    saving_dir: string path to directory to save figures under. Default=./
    show_all: bool whether or not to show all graphs. Default=False
    compress: bool whether or not to compress output files. If True, all subsequnce plots are compressed. Default=True
    hide_hybrids: bool whether or not to hide hybrid proteins in plots. Default=True
RETURNS: 
    None
'''
def plot_experiment(exp, agg_func='sum', show_all=False, saving_dir='./', compress=True, hide_hybrids=True):
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
    print('Generating peptide score plots...')
    plot_peptide_scores(exp, agg_func=agg_func, saving_dir=saving_dir, show_all=show_all, compress=compress, hide_hybrids=hide_hybrids)
    print('Finished.')

    # Plot the ranking of the corect 
    print('Generating score ranking plots...')
    peptide_plotting.plot_score_rankings(exp, save_dir=saving_dir + 'ranking_plots/', show_all=show_all)
    print('Finished.')

    # Plot score distributions vs protein sequence
    print('Generating score distributions vs protein sequences...')
    plot_protein_summary(exp, saving_dir=saving_dir, show_all=show_all, compress=compress)
    print('Finished.')