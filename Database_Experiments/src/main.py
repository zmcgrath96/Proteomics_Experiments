import os
import argparse
import json
from time import time
import database
from spectra import gen_spectra_files
from scoring import score_peptides
from sequences import peptides, proteins
from analysis import plotting, analyze_experiment
from utils import __file_exists, __make_valid_dir_string, __make_dir, __is_json

''' old hybrid "sequence": "ALYLVCGELYTSRV", 
    second hybid: GFFYTPKEANIR
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
    input_file = args.input_file
    num_peptides = args.num_peptides
    num_hybs = args.num_hybrids
    agg_func = args.agg_func
    show_all = args.show_all
    save_dir = __make_valid_dir_string(args.save_dir)
    __make_dir(save_dir)
    min_length = args.min_length
    max_length = args.max_length
    mix = args.mix
    compress = args.compress
    hide_hybs = args.hide_hybs
    digest = args.digest

    start_time = time()
    
    if not __is_json(input_file):
        # load in a list of proteins from a source file
        prots = proteins.load_proteins(input_file) 

        # make hybrid proteins
        hyb_prots = proteins.generate_hybrids(prots, num_hybs, min_contribution=max_length)
    
        # create peptides
        non_hybrid_peps = peptides.gen_peptides(prots, num_peptides, min_length=min_length, max_length=max_length, digest=digest)

        # create hybrid peptides
        hyb_peps = peptides.gen_peptides(hyb_prots, num_hybs, min_length=min_length, max_length=max_length, digest=digest, hybrid_list=True)

        # combine them for later use
        print('Combining hybrid and non hybrid lists...')
        all_proteins_raw = prots + hyb_prots
        all_proteins_cleaned = prots + [{'name': x['name'] , 'sequence': x['protein']} for x in hyb_prots]
        all_peptides_raw = non_hybrid_peps + hyb_peps
        all_peptides_cleaned = [{'name': x['peptide_name'], 'sequence': x['peptide_sequence']} for x in non_hybrid_peps] + [{'name': x['peptide_name'], 'sequence': x['peptide_sequence']} for x in hyb_peps]
        print('Done')

        # create database files
        print('Generating fasta databases...')
        fasta_databases = database.generate(all_peptides_cleaned, save_dir=save_dir)
        print('Done')

        # create spectrum files
        print('Generating spectra files...')
        spectra_files = gen_spectra_files.generate(all_proteins_cleaned, defaults['window_sizes'], save_dir=save_dir, compress=compress)
        print('Done.')

        # run scoring algorithm on database and k-mers
        print('Scoring...')
        score_output_files = score_peptides.score_peptides(spectra_files, fasta_databases, save_dir, compress=compress, crux_search=False, path_to_crux_cmd=defaults['crux_cmd'])
        print('Done.')

        # save scores to json
        print('Analyzing Experiment...')
        analyze_experiment.analyze(all_proteins_raw, all_peptides_raw, score_output_files, {**vars(args), **defaults}, predicting_agg_func=agg_func, saving_dir=save_dir, mix_in_hybrids=mix, show_all=show_all, compress=compress, hide_hybrids=hide_hybs)
        print('Done.')

    else:
        print('Loading experiment file...')
        experiment_json = json.load(open(input_file, 'r'))
        print('Finished loading experiment')
        plotting.plot_experiment(experiment_json, agg_func='', show_all=show_all, saving_dir=save_dir, compress=compress, hide_hybrids=hide_hybs)


    print('Finished experiment. Time to complete: {} seconds'.format(time() - start_time))

    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Entry file for the database experiments')
    parser.add_argument('input_file', type=str, metavar='F', help='Path to a file (csv or fasta) with protein name and sequences OR path to an experiment json file (with .json extension) to load and generate plots')
    parser.add_argument('--num-peptides', dest='num_peptides', type=int, default=50, help='Number of peptides to generate as the fake sample. Default=50')
    parser.add_argument('--num-hybrids', dest='num_hybrids', type=int, default=10, help='Number of hybrid proteins and peptides to generate. Default=10')
    parser.add_argument('--aggregate-function', dest='agg_func', type=str, default='sum', help='Which aggregation function to use for combining k-mer scores. Pick either sum or product. Default=sum')
    parser.add_argument('--show-all-graphs', dest='show_all', type=bool, default=False, help='Show all the graphs generated. Will save to directory either way. Default=False.')
    parser.add_argument('--output-dir', dest='save_dir', type=str, default='./', help='Directory to save all figures. Default=./')
    parser.add_argument('--min-length', dest='min_length', type=int, default=3, help='Minimum length peptide to create. Default=3')
    parser.add_argument('--max-length', dest='max_length', type=int, default=20, help='Maximum length peptide to create. Cuts from N terminus (left) side. Default=20')
    parser.add_argument('--measure-func', dest='m_func', type=str, default='average', help='Measuring function for determining the top n proteins. Options are: sum, average, max. Default=average')
    parser.add_argument('--peptide-file', dest='d_file', type=str, default='', help='Peptides from a past experiment. Default=None')
    parser.add_argument('--mix-prots', dest='mix', type=bool, default=True, help='Whether or not to also use hybrid proteins when calculating scores. Default=True')
    parser.add_argument('--compress', dest='compress', type=bool, default=True, help='Compress spectra files while generating them. Default=True')
    parser.add_argument('--hide-hybrid-prots', dest='hide_hybs', type=bool, default=False, help='When plotting hybrid peptides, hide the hybrid protein and only show results from normal proteins. Default=True')
    parser.add_argument('--digest', dest='digest', type=str, default='random', help='Type of digest to perform. Options are <random, trypsin>. Default=random')
    args = parser.parse_args()
    main(args)
    