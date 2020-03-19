import os
import sys
import argparse
import json
from time import time
from spectra import gen_spectra_files
from scoring import search_experiment
from sequences import peptides, proteins
from analysis import experiment
from plotting import plotting
from utils import file_exists, make_valid_dir_string, make_dir, is_json, is_fasta
from summarize import summary

''' 
    IMPORT DEFAULTS
'''
cwd = '/'.join(str(os.path.dirname(os.path.realpath(__file__))).split('/')[:-1])
default_json_file = cwd + '/defaults.json'
defaults = None
with open(default_json_file, 'r') as o:
    defaults = json.load(o)
'''
    END IMPORT DEFAULTS
'''
'''
    CONSTANTS
'''
starting_positions = ['b', 'p', 'a', 'sc', 'su']
'''
    END CONSTANTS
'''

#########################################################
#               MAIN HELPER FUNCTIONS
#########################################################
def str2bool(v):
    if isinstance(v, bool):
       return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

def is_correct_file(pos, file):
    if pos == 'b':
        return is_fasta(file), '.fasta'
    else: 
        return is_json(file), '.json'

def clean_prots(prots):
    cleaned = []
    for p in prots:
        if 'protein' in p:
            cleaned.append({
                'name': p['name'],
                'sequence': p['protein']
            })
        else:
            cleaned.append(p)
    return cleaned

def clean_peps(peps):
    return [{'name': p['peptide_name'], 'sequence': p['peptide_sequence']} for p in peps]

#########################################################
#           END MAIN HELPER FUNCTIONS
#########################################################

def main(args):
    '''
    STEPS
    1. Generate hybrid proteins and sequences (or load from file)
    2. Load proteins and generate (or load) peptides from digestion of proteins
    3. Save all peptides as .FASTA files
    4. For each protein, create k-mers and save to mzml files
    5. Score all the k-mers against all the peptides
    6. Perform analysis on these scores
    '''
    # input files
    input_file = args.input_file
    # output parameters
    save_dir =make_valid_dir_string(args.save_dir)
    make_dir(save_dir)
    # peptide/protein parameters
    num_peptides = args.num_peptides
    num_hybs = args.num_hybrids
    min_length = args.min_length
    max_length = args.max_length
    digest = args.digest
    len_dist = args.len_dist
    # analysis/plotting parameters
    agg_func = args.agg_func
    score_func = args.score_func
    show_all = args.show_all
    plot_pep_scores = args.plot_pep_scores
    plot_pep_ranks_length = args.plot_pep_ranks_length
    plot_pep_ranks_prot = args.plot_pep_ranks_prot
    # performance parameters
    compress = args.compress
    # flow control parameters
    start_pos = args.start_pos.lower() 
    start_pos = start_pos if start_pos in starting_positions else 'b'

    # check to see if the correct file was given
    c, ftype = is_correct_file(start_pos, input_file)
    if not c:
        print('Wrong input file. For starting position "{}" a "{}" is needed'.format(start_pos, ftype))
        sys.exit()

    start_time = time()
    experiment_json_file = None
    exp_json = None
    
    if start_pos == 'b':
        # load in a list of proteins from a source file
        prots = proteins.load_proteins(input_file) 

        # make hybrid proteins
        hyb_prots = proteins.generate_hybrids(prots, num_hybs, min_contribution=max_length)
    
        # create peptides
        non_hybrid_peps = peptides.gen_peptides(prots, num_peptides, min_length=min_length, max_length=max_length, digest=digest, dist=len_dist)

        # create hybrid peptides
        hyb_peps = peptides.gen_peptides(hyb_prots, num_hybs, min_length=min_length, max_length=max_length, digest=digest, hybrid_list=True)

        all_proteins_raw = prots + hyb_prots
        all_peptides_raw = non_hybrid_peps + hyb_peps

        print('Saving experiment...')
        experiment_json_file = experiment.save_experiment(all_proteins_raw, all_peptides_raw, {**vars(args), **defaults}, saving_dir=save_dir)
        print('Done.')

    if start_pos == 'sc' or start_pos == 'b':

        experiment_json_file = experiment_json_file if experiment_json_file is not None else input_file
        start_pos == 's' and print('Loading experiment file...')
        exp_json = json.load(open(experiment_json_file, 'r'))
        start_pos == 's' and print('Finished loading experiment')

        # combine them for later use
        print('Combining hybrid and non hybrid lists...')
        all_proteins_raw = exp_json['experiment_info']['proteins']
        all_proteins_cleaned = clean_prots(all_proteins_raw)
        all_peptides_raw = exp_json['experiment_info']['peptides']
        all_peptides_cleaned = clean_peps(all_peptides_raw)
        print('Done')

        # create database files
        print('Generating fasta databases...')
        fasta_databases = proteins.generate_databases(all_peptides_cleaned, save_dir=save_dir)
        print('Done')

        # create spectrum files
        print('Generating spectra files...')
        spectra_files = gen_spectra_files.generate(all_proteins_cleaned, defaults['window_sizes'], save_dir=save_dir, compress=compress)
        print('\nDone.')

        # run scoring algorithm on database and k-mers
        print('Scoring...')
        score_output_files = search_experiment.score_peptides(spectra_files, fasta_databases, save_dir, compress=compress, score_func=score_func, path_to_crux_cmd=defaults['crux_cmd'])
        print('\nDone.')

        # save the experiment in a json file
        print('Saving scores...')
        experiment_json_file = experiment.save_experiment(all_proteins_raw, all_peptides_raw, {**vars(args), **defaults}, files=score_output_files, saving_dir=save_dir)
        print('Done.')

    if start_pos == 'a' or start_pos == 'b' or start_pos == 'sc':
        # load experiment file
        print('Loading experiment...')
        experiment_json_file = experiment_json_file if experiment_json_file is not None else input_file
        exp_json = json.load(open(experiment_json_file, 'r'))
        print('Finished loading experiment file.')

        # Perform aggregations and k-mer ranksings
        print('Analyzing Experiment...')
        experiment_json_file, exp_json = experiment.analyze(exp_json, predicting_agg_func=agg_func, saving_dir=save_dir)
        print('\nDone.')

    if start_pos == 'a' or start_pos == 'b' or start_pos == 'sc' or start_pos == 'p':
        # check to see if exp has been loaded
        experiment_json_file = input_file if experiment_json_file is None else experiment_json_file
        if exp_json is None:
            print('Loading experiment file...')
            exp_json = json.load(open(experiment_json_file, 'r'))
            print('Finished loading experiment')
        # plot experiment
        plotting.plot_experiment(exp_json, agg_func='', show_all=show_all, saving_dir=save_dir, compress=compress, plot_pep_scores=plot_pep_scores, plot_pep_ranks_len=plot_pep_ranks_length, plot_pep_ranks_prot=plot_pep_ranks_prot)

    # check to see if exp has been loaded
    experiment_json_file = input_file if experiment_json_file is None else experiment_json_file
    if exp_json is None:
        print('Loading experiment file...')
        exp_json = json.load(open(experiment_json_file, 'r'))
        print('Finished loading experiment')
    print('Generating Sumamry...')
    summary.make_summary(exp_json, output_dir=save_dir)
    print('Done.')

    print('Finished experiment. Time to complete: {} seconds'.format(time() - start_time))

    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Entry file for the database experiments')
    parser.add_argument('input_file', type=str, metavar='F', help='Path to a fasta with protein name and sequences OR path to an experiment json file (with .json extension) to load')
    parser.add_argument('--start', type=str, dest='start_pos', default='b', help='Program position to start at. Four options: \nb --- beggining. Start from beginning with a .fasta file with proteins.\nsc --- scoring. Use proteins from an older generation and score them.\na --- analysis. Start at the analysis step with an old experiment json file.\np --- plotting. Use an old analysis from an experiment json and generate new plots.\nsu --- summary. Generate a summay file for an experiment.\n Default=b')
    parser.add_argument('--num-peptides', dest='num_peptides', type=int, default=50, help='Number of peptides to generate as the fake sample. Default=50')
    parser.add_argument('--num-hybrids', dest='num_hybrids', type=int, default=2, help='Number of hybrid proteins to generate. Default=2')
    parser.add_argument('--aggregate-function', dest='agg_func', type=str, default='sum', help='Which aggregation function to use for combining k-mer scores. Pick either sum or product. Default=sum')
    parser.add_argument('--show-all-graphs', dest='show_all', type=str2bool, default=False, help='Show all the graphs generated. Will save to directory either way. Default=False.')
    parser.add_argument('--output-dir', dest='save_dir', type=str, default='~/', help='Directory to save all figures. Default=~/')
    parser.add_argument('--min-length', dest='min_length', type=int, default=4, help='Minimum length peptide to create. Default=4')
    parser.add_argument('--max-length', dest='max_length', type=int, default=35, help='Maximum length peptide to create. Cuts from N terminus (left) side. Default=35')
    parser.add_argument('--length-dist', dest='len_dist', type=str, default='beta', help='Peptide length distribution to create. Options are: beta, random. Defualt=beta (based on expeimental).')
    parser.add_argument('--measure-func', dest='m_func', type=str, default='average', help='Measuring function for determining the top n proteins. Options are: sum, average, max. Default=average')
    parser.add_argument('--compress', dest='compress', type=str2bool, default=True, help='Compress spectra files while generating them. Default=True')
    parser.add_argument('--digest', dest='digest', type=str, default='random', help='Type of digest to perform. Options are <random, trypsin>. Default=random')
    parser.add_argument('--score-function', dest='score_func', type=str, default='custom', help='Type of scoring function to use. Options are "custom" or "crux". Default=custom')
    parser.add_argument('--plot-peptide-scores', dest='plot_pep_scores', type=str2bool, default=False, help='Determines whether to generate plots for every peptide score. Default=False')
    parser.add_argument('--plot-peptide-rankings-length', dest='plot_pep_ranks_length', type=str2bool, default=False, help='Determines whether to generate plots for peptide rankings distributions by length. Default=False')
    parser.add_argument('--plot-peptide-rankings-protein', dest='plot_pep_ranks_prot', type=str2bool, default=True, help='Determines whether to generate plots for peptide scoring ranks agaisnt proteins. Default=True')
    args = parser.parse_args()
    main(args)
    