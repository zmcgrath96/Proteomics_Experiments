import os
import argparse
import json
from database import gen_db
from spectra import gen_spectra_files
from scoring import score_peptides
from plotting import sequence_and_score
from sequences import peptides

''' old hybrid "sequence": "ALYLVCGELYTSRV", 
    second hybid: GFFYTPKEANIR
    IMPORT DEFAULTS
'''
cwd = os.path.dirname(os.path.realpath(__file__))
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
    1. Perform digest to get peptides
    2. Save all peptides (generated and hybrid) as FASTA files -- return the file names/paths
    3. Look at all proteins and generate k-mers and save as mzML files -- return the file names/paths
    4. For each peptide:
        For each protein: 
            For each k-mer length:
                score k-mers against peptide
            keep an aggregation of all k-mer scores 
        Compare how aggregations of each protein k-mer does, keep good ones (ones that aren't noise)
    5. Plot relevant scores of proteins
    '''
    experiment = 'fractionated' if 'fractionated' in str(args.experiment).lower() else 'flipped'
    num_peptides = args.num_peptides
    agg_func = args.agg_func
    show_all = args.show_all
    save_dir = args.save_dir
    '''
        SETUP ARGUMENTS FOR EACH STEP
    '''
    # load sequences once instead of all the time
    sequences_json = cwd + '/sequences.json'
    sequences = None
    with open(sequences_json, 'r') as seqfile:
        sequences = json.load(seqfile)
    
    db_args = {
        'experiment': experiment, 
        'path': cwd + '/' + defaults['save_dirs'][experiment], 
        'name': defaults['database_names'][experiment], 
        'window_sizes': defaults['window_sizes'], 
        'prefix': defaults['database_name_prefix'][experiment], 
        'sequences_dict': sequences, 
        'peptide_index': defaults['peptide_index']
        } 
    spectra_args = {
        'experiment': experiment, 
        'path': cwd + '/' + defaults['save_dirs'][experiment], 
        'name': defaults['spectra_names'][experiment], 
        'window_sizes': defaults['window_sizes'], 
        'title_prefix': defaults['spectrum_title_prefix'][experiment],
        'sequences_dict': sequences,
        'peptide_index': defaults['peptide_index']
        } 
    '''
        END ARGUMENT SETUP
    '''
    # create peptides
    peptides.gen_peptides(sequences, num_peptides, peptide_index=defaults['peptide_index'])
    # create database files
    fasta_databases = gen_db.generate(db_args)
    # create spectrum files
    spectra_files = gen_spectra_files.generate(spectra_args)
    # run scoring algorithm on database and k-mers
    score_output_files = score_peptides.score_peptides(spectra_files, fasta_databases, defaults['crux_cmd'], cwd + '/crux_output')
    # filter and plot scores
    protein_names = [x['name'] for x in sequences['sample']['proteins']]
    sequence_and_score.plot_experiment(experiment, score_output_files, protein_names, 'peptide', num_peptides, 'hybrid', agg_func=agg_func, show_all=show_all, saving_dir=save_dir)

    print('Finished.')
    print('===================================')
    print('SUMMARY\n\nRan {} experiment\nFASTA Databses:\n {}\n\nSpctra files: {}\n\nOutput files created: {}'.format(experiment, ', '.join(fasta_databases), ', '.join(spectra_files), ', '.join(score_output_files)))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Entry file for the database experiments')
    parser.add_argument('experiment', metavar='E', type=str, help='The experiment to run. Options are: \nflipped, fractionated\n. Defualts to flipped')
    parser.add_argument('--num-peptides', dest='num_peptides', type=int, default=49, help='Number of peptides to generate as the fake sample. Default 49')
    parser.add_argument('--aggregate-function', dest='agg_func', type=str, default='sum', help='Which aggregation function to use for combining k-mer scores. Pick either sum or product. Default sum')
    parser.add_argument('--show-all-graphs', dest='show_all', type=bool, default=False, help='Show all the graphs generated. Defaults to False. Will save to directory either way.')
    parser.add_argument('--save-directory', dest='save_dir', type=str, default='./', help='Directory to save all figures. Default is ./')
    args = parser.parse_args()
    main(args)
    