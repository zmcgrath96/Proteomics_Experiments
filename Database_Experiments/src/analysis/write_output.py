import json
from utils import __make_valid_text_file, __make_valid_json_file
from analysis import analysis_utils

            
'''write_raw_json

DESC:
    save the score data in a json file in dictionary or list form
PARAMS:
    scores: a dictionary of scores to save
    file: string a file path to save all the data in 
RETURNS:
    none
'''
def write_raw_json(file, scores):
    file = __make_valid_json_file(file)
    if not type(scores) == dict and not type(scores) == list:
        raise Exception('scores should be type dict or list. is {}'.format(type(scores)))
    with open(file, 'w') as o:
        json.dump(scores, o)