from os import walk
# get the score files from 
path_to_score_files = '/Users/zacharymcgrath/Desktop/Experiment output/NOD2_E3/search_output/'

score_files = []

for (dirpath, dirnames, filenames) in walk(path_to_score_files):
    for fname in filenames:
        score_files.append(path_to_score_files + fname)
    break


import json
print('loading json')
experiment_json_file = '/Users/zacharymcgrath/Desktop/Experiment output/NOD2_E3/experiment_data.json'
exp = json.load(open(experiment_json_file, 'r'))
print('done loading json')

prots = exp['experiment_info']['proteins']
peps = exp['experiment_info']['peptides']

import sys
sys.path.append('/Users/zacharymcgrath/Documents/Layer_Research/Proteomics_Experiments/Database_Experiments/src')

from analysis import experiment

print('saving experiment')
experiment.save_experiment(prots, peps, {}, files=score_files)