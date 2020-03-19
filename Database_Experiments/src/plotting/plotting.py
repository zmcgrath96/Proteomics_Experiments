from matplotlib import pyplot as plt 
import json
import itertools
import os
import numpy as np
from copy import deepcopy
from utils import __get_related_files, __make_dir, __make_valid_dir_string, __gzip_dir, __split_exp_by_ion, experiment_has_ion_types
from analysis.analysis_utils import get_top_n, __get_argmax_max
from analysis.score_utils import align_scan_pos, get_scores_scan_pos_label
from plotting import peptide_plotting, protein_plotting

####################################################
#               CONSTANTS
####################################################
agg_funcs = ['product', 'sum', 'z_score_sum']

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
Inputs:
    scores: dict of list of floats from the analysis
kwargs:
    agg_func: str name of the aggregation function. If none is given, looks for known aggregation functions. Default=''
Outputs:
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
            return pep['starting_position']
    return -1

'''__aggregate_ranks

DESC:
    sort out the ranks for each peptide from the experiment json
    These ranks are the ranks of the aggregate function 
    from the correct protein at the correct position
Inputs:
    exp: dictionary that holds the experiment information
Outputs:
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
Inputs:
    exp: dict the json used to save all information 
kwargs:
    saving_dir: string path to directory to save figures under. Default=./
    show_all: bool whether or not to show all graphs. Default=False
    compress: bool whether or not to compress output files. If True, all subsequnce plots are compressed. Default=True
Outputs:
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
Inputs:
    exp: dict the json used to save all information
kwargs:
    agg_func: string aggregate function to use. Default=sum
    saving_dir: string path to directory to save figures under. Default=./
    show_all: bool whether or not to show all graphs. Default=False
    compress: bool whether or not to compress output files. If True, all subsequnce plots are compressed. Default=True
Outputs: 
    None
'''
def plot_peptide_scores(exp, agg_func='sum', saving_dir='./', show_all=False, compress=True):
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
        peptide_plotting.plot_subsequence(agg_scores, title=str(peptide), save_dir=pep_saving_dir, show_graph=show_all, agg_func=agg_func, sequence_info=info, compress=compress)
    
def plot_generator(exp: dict, agg_func='sum', show_all=False, saving_dir='./', compress=True, plot_pep_scores=True, plot_pep_ranks_len=True, plot_pep_ranks_prot=True) -> None:
    '''
    Generate various plots for the experiment

    Inputs:
        exp: dict the json used to save all information
    kwargs:
        agg_func: string aggregate function to use. Default=sum
        saving_dir: string path to directory to save figures under. Default=./
        show_all: bool whether or not to show all graphs. Default=False
        compress: bool whether or not to compress output files. If True, all subsequnce plots are compressed. Default=True
    Outputs: 
        None
    '''
     #create the saving directory
    saving_dir = __make_valid_dir_string(saving_dir)
    __make_dir(saving_dir)        

    # Plot the kmer scores and score aggregations
    if plot_pep_scores:
        print('Generating peptide score plots...')
        plot_peptide_scores(exp, agg_func=agg_func, saving_dir=saving_dir, show_all=show_all, compress=compress)
        print('Finished.')

    # Plot the ranking of the corect 
    if plot_pep_ranks_len:
        print('Generating score ranking plots...')
        peptide_plotting.plot_score_rankings(exp, save_dir=saving_dir + 'ranking_plots/', show_all=show_all)
        print('Finished.')

    # Plot score distributions vs protein sequence
    if plot_pep_ranks_prot:
        print('Generating score distributions vs protein sequences...')
        plot_protein_summary(exp, saving_dir=saving_dir, show_all=show_all, compress=compress)
        print('Finished.')

def plot_experiment(exp: dict, agg_func='sum', show_all=False, saving_dir='./', compress=True, plot_pep_scores=True, plot_pep_ranks_len=True, plot_pep_ranks_prot=True) -> None:
    '''
    Generate various plots for the experiment
    
    Inputs:
        exp: dict the json used to save all information
    kwargs:
        agg_func: string aggregate function to use. Default=sum
        saving_dir: string path to directory to save figures under. Default=./
        show_all: bool whether or not to show all graphs. Default=False
        compress: bool whether or not to compress output files. If True, all subsequnce plots are compressed. Default=True
    Outputs: 
        None
    '''
    print('Generating plots...')
    if experiment_has_ion_types(exp):
        for ion in ['b', 'y']:
            ion_exp = __split_exp_by_ion(exp, ion)
            plot_generator(
                ion_exp, 
                agg_func=agg_func, 
                show_all=show_all, 
                saving_dir=__make_valid_dir_string(saving_dir) + 'ion_{}'.format(ion), 
                compress=compress, 
                plot_pep_scores=plot_pep_scores, 
                plot_pep_ranks_len=plot_pep_ranks_len,
                plot_pep_ranks_prot=plot_pep_ranks_prot
            )
    else:
        plot_generator(
            exp, 
            agg_func=agg_func, 
            show_all=show_all, 
            saving_dir=saving_dir, 
            compress=compress, 
            plot_pep_scores=plot_pep_scores, 
            plot_pep_ranks_len=plot_pep_ranks_len,
            plot_pep_ranks_prot=plot_pep_ranks_prot
        )
    print('Finished generating all plots.')