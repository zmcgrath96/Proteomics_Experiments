import argparse
import pandas as pd 
import json
from copy import deepcopy
import aggregations

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

def __get_pep_row_by_seq(df, seq, col):
    return df.loc[df[col] == seq]

def __get_pep_start(pep_seq):
    global digest
    pep_row = __get_pep_row_by_seq(digest, pep_seq, dig_pep_seq)
    return int(pep_row[dig_pep_start])

def __get_pep_name(seq):
    global digest, exp
    possible_names = [x['peptide_name'] for x in exp['header']['peptides']]
    this_name = str(digest.loc[digest[dig_pep_seq] == seq][dig_name_col])
    for name in possible_names:
        if name in this_name:
            return name

def __get_predicted_parents(pep_seq):
    global exp
    peps = exp['header']['peptides']
    name = None 
    for pep in peps:
        if pep['peptide_sequence'] == pep_seq:
            name = pep['peptide_name']
            break 
    parents = exp['experiment'][name]['predicted_parents']
    p_names = [x['protein_name'] for x in parents]
    return p_names

def __get_start_score(pep_seq, agg):
    global exp
    agg_f = aggregations.__z_score_sum if 'z_score_sum' in agg.lower() else (aggregations.__product if 'product' in agg.lower() else aggregations.__sum)
    start = __get_pep_start(pep_seq)
    name = __get_pep_name(pep_seq)
    prots = exp['experiment'][name]
    aggs = {}
    for prot in prots:
        if prot == 'predicted_parents':
            continue
        mers = [mer for _, mer in prots[prot].items()]
        aggs[prot] = agg_f(mers)

    max_s = -100
    max_prot = ''
    for prot, a in aggs.items():
        if start > len(a):
            continue
        if a[start] > max_s:
            max_s = a[start]
            max_prot = prot
    return max_s, max_prot

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
        row = __get_pep_row_by_seq(crux_output, pep, crux_seq_col)
        rel.append(row)
    return pd.concat(rel)

def __gen_comparison(relevant, agg_func):
    form = {
        'sequence': None, 
        'crux_predicted_parents': None, 
        'crux_score': None,
        'flipped_score_start': None,
        'flipped_prot_highest_start': None,
        'flipped_predicted_parents': None,
    }
    for _, row in relevant.iterrows():
        t = deepcopy(form)
        t['sequence'] = row[crux_seq_col]
        t['crux_predicted_parents'] = row[crux_parent_col]
        t['crux_score'] = row[crux_score]
        start, max_prot_at_start = __get_start_score(t['sequence'], agg_func) 
        t['flipped_score_start'] = start
        t['flipped_prot_highest_start'] = max_prot_at_start
        t['flipped_predicted_parents'] = __get_predicted_parents(t['sequence'])
        print(row['c1'], row['c2'])


####################################################
#            end utilities
####################################################

def compare(args):
    __load_files(args)
    relevant = __get_relevant_rows()
    __gen_comparison(relevant, args.agg_func)

    
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('crux_output', type=str, metavar='CO', help='The output file path from the traditional search method')
    parser.add_argument('digestion', type=str, metavar='D', help='The digestion TSV file path from the flipped experiment')
    parser.add_argument('experiment_json', type=str, metavar='EJ', help='The json file path generated from the experiment')
    parser.add_argument('agg_func', type=str, metavar='AF', help='The aggregate function used in the experiment')
    args = parser.parse_args()
    compare(args)