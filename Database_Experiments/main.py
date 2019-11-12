import os
import argparse
import json
from database import gen_db, load_dbs
from spectra import gen_spectra_files, load_spectra

'''
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
    1. Generate databases and get the file names
    2. Generate the spectra and get the file names
    3. Load the database into memory
    4. Load the spectra into memory
    5. Score each spectra against the database
    6. Plot a graph of the score against the sequence position
    '''
    experiment = 'fractionated' if 'fractionated' in str(args.experiment).lower() else 'flipped'
    '''
        SETUP ARGUMENTS FOR EACH STEP
    '''
    db_args = {
        'experiment': experiment, 
        'path': cwd + '/' + defaults['save_dirs'][experiment], 
        'name': defaults['database_names'][experiment], 
        'window_sizes': defaults['window_sizes'], 
        'prefix': defaults['database_name_prefix'][experiment], 
        'sequences_json': cwd + '/sequences.json'
        } 
    spectra_args = {
        'experiment': experiment, 
        'path': cwd + '/' + defaults['save_dirs'][experiment], 
        'name': defaults['spectra_names'][experiment], 
        'window_sizes': defaults['window_sizes'], 
        'title_prefix': defaults['spectrum_title_prefix'][experiment],
        'sequences_json': cwd + '/sequences.json'
        } 
    '''
        END ARGUMENT SETUP
    '''
    fasta_databases = gen_db.generate(db_args)
    spectra_files = gen_spectra_files.generate(spectra_args)
    # reads are done via the pyteomics module
    loaded_databases = load_dbs.load_dbs(fasta_databases)
    loaded_spectra = load_spectra.load_spectra(spectra_files)



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Entry file for the database experiments')
    parser.add_argument('experiment', metavar='E', type=str, help='The experiment to run. Options are: \nflipped, fractionated\n. Defualts to flipped')
    args = parser.parse_args()
    main(args)
    