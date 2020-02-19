from copy import deepcopy
import json
import utils
from analysis.write_output import write_raw_json
from analysis import score_utils
from analysis.aggregations import __z_score_sum, __sum, __product
from analysis.analysis_utils import get_top_n_prots
from analysis.plotting import plot_experiment

#######################################################
#                   CONSTANTS
#######################################################
experiment_json_file_name = 'experiment_data.json'

EXPERIMENT_ENTRY = 'experiment'
EXPERIMENT_HEADER = 'experiment_info'
EXPERIMENT_PROTEIN_HEADER = 'proteins'
EXPERIMENT_PEPTIDE_HEADER = 'peptides'
EXPERIMENT_ARGUMENT_HEADER = 'arguments'
EXPERIMENT_PARENT_PREDICTION = 'predicted_parents'
SAMPLE_ENTRY = 'sample'
SAMPLE_PROTEINS = 'proteins'
SAMPLE_PROTEIN_ANALYSIS = 'analysis'

HYBRID_SEACH_STRING = 'HYBRID'
#####################################################
#               END CONSTANTS
#####################################################

#####################################################
#               GLOBALS
#####################################################
experiment_json = {
    EXPERIMENT_HEADER: {
        EXPERIMENT_PROTEIN_HEADER: None, 
        EXPERIMENT_PEPTIDE_HEADER: None,
        EXPERIMENT_ARGUMENT_HEADER: None
    }, 
    EXPERIMENT_ENTRY: {

    }
}

#####################################################
#               END GLOBALS
#####################################################

#####################################################
#                "PRIVATE" FUNCTIONS
#####################################################

'''
'''
def __get_k_number(file_name):
    left_side = str(file_name).split('vs')[0]
    no_u = left_side.split('_')
    k = int(no_u[-1]) if no_u[-1] is not '' else int(no_u[-2])
    return k

'''__add_header_info
DESC:
    add all header info to the experiment json
PARAMS:
    proteins: list of dictionaries of the form {name: str, sequence: str}
    peptides: list of dictionaries of the form 
    {    
        'peptide_name': str,
        'peptide_sequence': str,
        'parent_name': str,
        'parent_sequence': str,
        'start_index': int, 
        'end_index': int
    }
    args: dictionary parameters used when running the experiment
    json: dictionary object in which to save the header info
'''
def __add_header_info(proteins, peptides, args, json):
    json[EXPERIMENT_HEADER][EXPERIMENT_PROTEIN_HEADER] = deepcopy(proteins)
    json[EXPERIMENT_HEADER][EXPERIMENT_PEPTIDE_HEADER] = deepcopy(peptides)
    json[EXPERIMENT_HEADER][EXPERIMENT_ARGUMENT_HEADER] = deepcopy(args)

'''__add_subsequence_agg

DESC:
    adds the aggregation information to each subsequence part in the dictionary
PARAMS:
    subsequence_name: string name of the subsequence
    protein_names: names of all the proteins 
    json: dictionary object where all items are save
OPTIONAL:
    predicting_agg_func: str name of the function used for aggregation. Default=sum
RETURNS: 
    None
'''
def __add_subsequnce_agg(peptide_dict, predicting_agg_func='sum'):
    agg_func = __z_score_sum if 'z_score_sum' in predicting_agg_func.lower() else ( __product if 'product' in predicting_agg_func.lower() else __sum)

    subsequence_aggs = {}
    for prot_name, prot_scores in peptide_dict.items():
    
        if SAMPLE_PROTEIN_ANALYSIS == prot_name:
            continue
            
        agged = agg_func(prot_scores)
        peptide_dict[prot_name][predicting_agg_func] = deepcopy(agged)
        subsequence_aggs[prot_name] = deepcopy(agged)
    
    if SAMPLE_PROTEIN_ANALYSIS not in peptide_dict: 
        peptide_dict[SAMPLE_PROTEIN_ANALYSIS] = {}

    peptide_dict[SAMPLE_PROTEIN_ANALYSIS][EXPERIMENT_PARENT_PREDICTION] = get_top_n_prots(subsequence_aggs)
    return peptide_dict   

'''__find_kmer_rank

DESC:
    find how well the correct k-mer scores against all other k-mers of the same k
    1 based scores
PARAMS:
    correct_prot: str name of the correct protein
    peptide_analysis: dictionary containg all the peptide stuff from analysis
RETURNS:
    dictionary of ranks. entry is name of kmer or aggregate and rank is its inner kmer rank
'''
def __find_kmer_rank(correct_prot, correct_position, peptide_analysis):
    tagged_scores = {}
    ranking = {}
    for prot, kmers in peptide_analysis.items():
        if prot == SAMPLE_PROTEIN_ANALYSIS:
            continue

        # create tuples of the form (parent protein name, score, order number)
        ranking[prot] = {}
        for k, scores in kmers.items():
            if k not in tagged_scores:
                tagged_scores[k] = []
            ranking[prot][k] = None
            tagged = [(prot, score, i) for i, score in enumerate(scores)]
            tagged_scores[k] += tagged

    for k, tg in tagged_scores.items():
        tg.sort(reverse=True, key=lambda x: x[1])
        # if scores are all the same, don't increment the rank
        rank = 0
        last_score = None
        for score in tg:
            fscore = float(score[1])
            pos = int(score[2])
            if fscore != last_score:
                rank += 1
                last_score = fscore
            if pos == correct_position:
                ranking[score[0]][k] = rank
    
    return ranking[correct_prot]

'''__rank_pep

DESC:
    run through all the proteins against this peptide and rank the correct score at the right position
PARAMS: 
    json: dictionary object to add analysis to
    peptide: dictionary of all the peptide information
RETURNS:
    None
'''
def __rank_pep(json, peptide):
    ranking_dict = {}
    ranking_dict['correct_protein'] = peptide['parent_name']
    ranking_dict['ranks'] = __find_kmer_rank(peptide['parent_name'], peptide['start_index'], deepcopy(json[EXPERIMENT_ENTRY][peptide['peptide_name']]))
    ranking_dict['sequence'] = peptide['peptide_sequence']
    ranking_dict['sequence_length'] = len(peptide['peptide_sequence'])
    json[EXPERIMENT_ENTRY][peptide['peptide_name']][SAMPLE_PROTEIN_ANALYSIS]['ranks'] = ranking_dict

#####################################################
#               END "PRIVATE" FUNCTIONS
#####################################################

'''save_experiment

DESC:
    builds an initial dictionary to be saved that has basic peptide, protein and score info
PARAMS:
    PARAMS:
    proteins: list of dictionaries of the form {name: str, sequence: str}
    peptides: list of dictionaries of the form 
    {    
        'peptide_name': str,
        'peptide_sequence': str,
        'parent_name': str,
        'parent_sequence': str,
        'start_index': int, 
        'end_index': int
    }
    args: dictionary parameters used when running the experiment
OPTIONAL:
    files: list of str output files from the scoring algorithm
    saving_dir: str path to where the experiment file should be saved. Default=./
RETURNS:
    str path to experiment json file
'''
def save_experiment(proteins, peptides, args, files=None, saving_dir='./'):
    global experiment_json_file_name, experiment_json
    saving_dir = utils.__make_valid_dir_string(saving_dir)
    utils.__make_dir(saving_dir)

    # add header information to the json
    __add_header_info(proteins, peptides, args, experiment_json)
    protein_names = [x['name'] for x in proteins]

    # if no scoring files exist, just save the proteins and peptides
    if files is None or (isinstance(files, list) and len(files) == 0):
        write_raw_json(saving_dir + experiment_json_file_name, experiment_json)
        return(saving_dir + experiment_json_file_name)

    # go through each peptide
    for pc, pep in enumerate(peptides):
        print('Progress: {}%\r'.format(int((float(pc)/float(len(peptides))) * 100)), end='')
        # get the peptide related files
        pep_related = utils.__get_related_files(files, pep['peptide_name'])
        subsequence_dict = {}
        # get the protein information for each peptide
        for prot_name in protein_names:
            subsequence_dict[prot_name] = {}
            prot_with_subseq = utils.__get_related_files(pep_related, str(prot_name).lower())
            if prot_with_subseq is None or len(prot_with_subseq) == 0: 
                print('No files scoring {} against {} were found. Skipping'.format(prot_name, pep['peptide_name']))

            this_prot_kmers = []
            for f in prot_with_subseq:
                scores, _, _ = score_utils.get_scores_scan_pos_label(f)
                k = 'k=' + str(__get_k_number(f))
                subsequence_dict[prot_name][k] = scores
                this_prot_kmers.append(scores)
                
        experiment_json[EXPERIMENT_ENTRY][pep['peptide_name']] = deepcopy(subsequence_dict)

    write_raw_json(saving_dir + experiment_json_file_name, experiment_json)
    return(saving_dir + experiment_json_file_name)
        

'''analyze

DESC:
    perform k-mer rankings and aggregations
PARAMS:
    exp: either a str to an experiment json file or a dictionary of a preloaded experiment file
OPTIONAL:
    predicting_agg_func: str name of the aggregation function to use. Default=sum
    saving_dir: str the name of the directory to save the experiment in. Default=./
RETURNS:
    str file path to the experiment json generated, dictionary of all experiment information
'''
def analyze(exp, predicting_agg_func='sum', saving_dir='./'):
    global experiment_json, experiment_json_file_name
    saving_dir = utils.__make_valid_dir_string(saving_dir)
    utils.__make_dir(saving_dir)

    # either load json from file or change its name
    if isinstance(exp, dict):
        experiment_json = exp
    elif isinstance(exp, str) and utils.__is_json(exp):
        experiment_json_file_name = exp 
        experiment_json = json.load(open(experiment_json_file_name, 'r'))
    
    # separate file name from full directory
    if len(experiment_json_file_name.split('/')) > 1:
        experiment_json_file_name = experiment_json_file_name.split('/')[-1]

    # perform analysis on all peptides
    peptides = experiment_json[EXPERIMENT_ENTRY]
    p_counter = 0
    peptide_header_list_idx = {}
    for i, d in enumerate(experiment_json[EXPERIMENT_HEADER][EXPERIMENT_PEPTIDE_HEADER]):
        peptide_header_list_idx[d['peptide_name']] = i

    for pep_name, peptide in peptides.items():
        p_counter += 1
        print('Analyzing peptide {}/{}[{}%]\r'.format(p_counter, len(peptides), int( (float(p_counter)/float(len(peptides))) *100 )), end='')
        peptide = __add_subsequnce_agg(peptide, predicting_agg_func=predicting_agg_func)
        peptide_info_dict = experiment_json[EXPERIMENT_HEADER][EXPERIMENT_PEPTIDE_HEADER][peptide_header_list_idx[pep_name]]
        __rank_pep(experiment_json, peptide_info_dict)
    
    # save to file
    write_raw_json(saving_dir + experiment_json_file_name, experiment_json)

    return(saving_dir + experiment_json_file_name), experiment_json
