import argparse
import pandas as pd 
import json
from copy import deepcopy

digest = None
crux_output = None 
exp = None
dig_name_col = 'peptide-name'
dig_pep_seq = 'peptide'
dig_pep_start = 'start-location'
crux_score = 'xcorr score'
crux_seq_col = 'sequence'
crux_parent_col = 'protein id'

####################################################
#               utilities
####################################################
def __load_json(f):
    j = None 
    with open(f, 'r') as i:
        j = json.load(i)
    return j

def __load_tsv(f):
    df = pd.read_csv(f, sep='\t')
    return df 

def __get_pep_row_by_seq(df, seq):
    return df.loc[df[dig_pep_seq] == seq]

def __get_pep_start(pep_seq):
    global digest
    pep_row = __get_pep_row_by_seq(digest, pep_seq)
    return int(pep_row[dig_pep_start])

def __load_files(args):
    global digest, crux_output, exp
    digest = __load_tsv(args.digestion)
    crux_output = __load_tsv(args.crux_output)
    exp = __load_json(args.experiment_json)

def __get_relevant_rows():
    global digest, crux_output
    in_sample = list(digest[dig_pep_seq])
    rel = []
    for pep in in_sample:
        rel.append(__get_pep_row_by_seq(crux_output, pep))
    return pd.concat(rel)

def __gen_comparison(relevant):
    form = {
        'sequence': None, 
        'crux_predicted_parents': None, 
        'crux_score': None,
        'flipped_score_start': None,
        'flipped_predicted_parents': None
    }
    for _, row in relevant.iterrows():
        t = deepcopy(form)
        t['sequence'] = row[crux_seq_col]
        t['crux_predicted_parents'] = row[crux_parent_col]
        t['crux_score'] = row[crux_score]
        t['flipped_score_start'] = __get_pep_start()
        print(row['c1'], row['c2'])


####################################################
#            end utilities
####################################################

def compare(args):
    __load_files(args)
    relevant = __get_relevant_rows()
    __gen_comparison(relevant)

    
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('crux_output', type=str, metavar='CO', help='The output file path from the traditional search method')
    parser.add_argument('digestion', type=str, metavar='D', help='The digestion TSV file path from the flipped experiment')
    parser.add_argument('experiment_json', type=str, metavar='EJ', help='The json file path generated from the experiment')
    parser.add_argument('agg_func', type=str, metavar='AF', help='The aggregate function used in the experiment')
    args = parser.parse_args()
    compare(args)