import os
import argparse
import json
from database import gen_db
from spectra import gen_spectra_files
from scoring import score_peptides
from sequence_generation import peptides
from analysis.plotting import plot_experiment
from analysis.analyze_experiment import analyze
from utils import __file_exists

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
    num_peptides = args.num_peptides
    agg_func = args.agg_func
    show_all = args.show_all
    save_dir = args.save_dir
    min_length = args.min_length
    max_length = args.max_length
    top_n = args.top_n
    n = args.n
    m_func = args.m_func
    test_seq = args.test_seq
    old_digest = args.d_file 
    '''
        SETUP ARGUMENTS FOR EACH STEP
    '''
    # load sequences once instead of all the time
    sequences_json = cwd + '/sequences.json' if not test_seq else cwd + '/test-sequences.json'
    sequences = None
    with open(sequences_json, 'r') as seqfile:
        sequences = json.load(seqfile)
    
    db_args = {
        'path': cwd + '/' + defaults['save_dirs'], 
        'name': defaults['database_names'], 
        'window_sizes': defaults['window_sizes'], 
        'prefix': defaults['database_name_prefix'], 
        'sequences_dict': sequences, 
        'peptide_index': defaults['peptide_index']
        } 
    spectra_args = {
        'path': cwd + '/' + defaults['save_dirs'], 
        'name': defaults['spectra_names'], 
        'window_sizes': defaults['window_sizes'], 
        'title_prefix': defaults['spectrum_title_prefix'],
        'sequences_dict': sequences,
        'peptide_index': defaults['peptide_index']
        } 
    '''
        END ARGUMENT SETUP
    '''
    # create peptides
    if not old_digest or old_digest == '' or not __file_exists(old_digest):
        peptides.gen_peptides(sequences, num_peptides, peptide_index=defaults['peptide_index'], min_length=min_length, max_length=max_length, save_dir=save_dir)
    else:
        peptides.load_peptides(sequences, old_digest, peptide_index=defaults['peptide_index']) 
    # create database files
    fasta_databases = gen_db.generate(db_args)
    # create spectrum files
    spectra_files = gen_spectra_files.generate(spectra_args)
    # run scoring algorithm on database and k-mers
    print('Scoring...')
    score_output_files = score_peptides.score_peptides(spectra_files, fasta_databases, defaults['crux_cmd'], save_dir + '/crux_output')
    print('Done.')
    # save scores to json
    protein_names = [x['name'] for x in sequences['sample']['proteins']]
    print('Saving experiment...')
    exp_json_path = analyze(score_output_files, protein_names, 'peptide', num_peptides, 'hybrid_db', sequences, saving_dir=save_dir, predicting_agg_func=agg_func, digestion_file=old_digest)
    print('Done.')
    # load the experiment and plot it
    plot_experiment(exp_json_path, agg_func=agg_func, show_all=show_all, saving_dir=save_dir, use_top_n=top_n, n=n, measure=m_func)

    print('Finished.')
    print('===================================')
    print('SUMMARY\n\nRan {} experiment\nFASTA Databses:\n {}\n\nSpctra files: {}\n\nOutput files created: {}'.format(experiment, ', '.join(fasta_databases), ', '.join(spectra_files), ', '.join(score_output_files)))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Entry file for the database experiments')
    parser.add_argument('--num-peptides', dest='num_peptides', type=int, default=49, help='Number of peptides to generate as the fake sample. Default=49')
    parser.add_argument('--aggregate-function', dest='agg_func', type=str, default='sum', help='Which aggregation function to use for combining k-mer scores. Pick either sum or product. Default=sum')
    parser.add_argument('--show-all-graphs', dest='show_all', type=bool, default=False, help='Show all the graphs generated. Will save to directory either way. Default=False.')
    parser.add_argument('--output-dir', dest='save_dir', type=str, default='./', help='Directory to save all figures. Default=./')
    parser.add_argument('--min-length', dest='min_length', type=int, default=3, help='Minimum length peptide to create. Default=3')
    parser.add_argument('--max-length', dest='max_length', type=int, default=20, help='Maximum length peptide to create. Cuts from N terminus (left) side. Default=20')
    parser.add_argument('--top-n', dest='top_n', type=bool, default=False, help='When recording how well a peptide scores against a protein, only use the top n proteins. Default=False')
    parser.add_argument('--n', dest='n', type=int, default=5, help='n to use if using --top-n. Default=5')
    parser.add_argument('--measure-func', dest='m_func', type=str, default='average', help='Measuring function for determining the top n proteins. Options are: sum, average, max. Default=average')
    parser.add_argument('--test-seq', dest='test_seq', type=bool, default=False, help='FOR TESTING ON SMALLER SEQUENCES. DEFAULT=False')
    parser.add_argument('--digestion-file', dest='d_file', type=str, default='', help='Digestion from a past experiment. Default=None')
    args = parser.parse_args()
    main(args)
    