import json
from utils import __make_valid_text_file, __make_valid_json_file
from data import analysis

SUMMARY_HEADER = '''
###############################################
                   SUMMARY 
###############################################
'''

TOP_5_HEADER = '''
TOP 5 SCORING PROTEINS
'''

RAW_HEADER = '''
RAW DATA
'''

'''write_summary

DESC: 
    save data from a plot to a text file
PARAMS:

'''
def write_summary(file, data, title=''):
    file = __make_valid_text_file(file)
    if type(data) is not dict: 
        raise Exception('scores should be type dict. is {}'.format(type(data)))
    with open(file, 'w') as o:
        o.write(str(title) + '\n')
        o.write(SUMMARY_HEADER)
        protein_line = 'proteins: {}\n'.format(', '.join(key for key in data))
        o.write(protein_line)
        top_5 = analysis.get_top_n(data)
        for i in range(len(top_5)):
            line = '{}:\n {},\n {}\n'.format(i, top_5[i][1], top_5[i][0])
            o.write(line)
    #     o.write(RAW_HEADER)
    # with open(file, 'a') as o:
    #     json.dump(data, file)
        
            
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