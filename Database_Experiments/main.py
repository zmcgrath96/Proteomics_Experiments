import os
import argparse
import json

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
    experiment = args.experiment
    db_args = {}
    spectra_args = {}
    if experiment == 'fractionated' or 'fractionated' in experiment:
        # do the fractionated experiments
        pass

    else:
        # do the flipped experiments 
        pass

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Entry file for the database experiments')
    parser.add_argument('experiment', metavar='E', type=str, help='The experiment to run. Options are: \nflipped, fractionated\n. Defualts to flipped')
    