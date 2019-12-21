from copy import deepcopy
import utils
from sequence_generation.digest import load_digest
from analysis.write_output import write_raw_json
from analysis import score_utils
from analysis.aggregations import __z_score_sum, __sum, __product
from analysis.analysis_utils import get_top_n_prots

#######################################################
#                   CONSTANTS
#######################################################
experiment_json_file_name = 'experiment_data.json'
EXPERIMENT_ENTRY = 'experiment'
EXPERIMENT_HEADER = 'header'
EXPERIMENT_PROTEIN_HEADER = 'proteins'
EXPERIMENT_PEPTIDE_HEADER = 'peptides'
EXPERIMENT_PARENT_PREDICTION = 'predicted_parents'
SAMPLE_ENTRY = 'sample'
SAMPLE_PROTEINS = 'proteins'
SAMPLE_HYBRID_ENTRY = 'hybrid'
SAMPLE_HYBRID_SEQUENCE = 'sequence'
SAMPLE_HYBRID_PARENT = 'full_parent'
SAMPLE_HYBRID_INDICES = 'parent_indices'
SAMPLE_PROTEIN_ANALYSIS = 'analysis'
#####################################################
#               END CONSTANTS
#####################################################

#####################################################
#               GLOBALS
#####################################################
experiment_json = {
    EXPERIMENT_HEADER: {
        EXPERIMENT_PROTEIN_HEADER: None, 
        EXPERIMENT_PEPTIDE_HEADER: []
    }, 
    EXPERIMENT_ENTRY: {

    }
}

digests = None
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
    json: dictionary object in which to save the header info
'''
def __add_header_info(proteins, peptides, json):
    json[EXPERIMENT_HEADER][EXPERIMENT_PROTEIN_HEADER] = deepcopy(proteins)
    json[EXPERIMENT_HEADER][EXPERIMENT_PEPTIDE_HEADER] = deepcopy(peptides)

'''__save_subsequence_info

DESC:
    gets and retrieves all information of the subsequence 
    Updates the global experiment_json object
PARAMS:
    subsequence_files: list of strings filepaths to files of the experiment
    subsequence_name: string name of the subsequence
    protein_names: names of all the proteins 
    json: dictionary object where all items are saved
RETURNS: 
    None
'''
def __save_subsequence_info(subsequence_files, subsequence_name, protein_names, json, predicting_agg_func='sum'):
    agg_func = __z_score_sum if 'z_score_sum' in predicting_agg_func.lower() else ( __product if 'product' in predicting_agg_func.lower() else __sum)
    subsequence_dict = {}
    subsequence_aggs = {}
    for prot_name in protein_names:
        subsequence_dict[prot_name] = {}
        prot_with_subseq = utils.__get_related_files(subsequence_files, str(prot_name).lower())
        if prot_with_subseq is None or len(prot_with_subseq) == 0: 
            print('No files scoring {} against {} were found. Skipping'.format(prot_name, subsequence_name))
        this_prot_kmers = []
        for f in prot_with_subseq:
            scores, _, _ = score_utils.__get_scores_scan_pos_label(f)
            k = 'k=' + str(__get_k_number(f))
            subsequence_dict[prot_name][k] = scores
            this_prot_kmers.append(scores)
        agged = agg_func(this_prot_kmers)
        subsequence_aggs[prot_name] = deepcopy(agged)
        subsequence_dict[prot_name][predicting_agg_func] = deepcopy(agged)
    
    if SAMPLE_PROTEIN_ANALYSIS not in subsequence_dict: 
        subsequence_dict[SAMPLE_PROTEIN_ANALYSIS] = {}
    subsequence_dict[SAMPLE_PROTEIN_ANALYSIS][EXPERIMENT_PARENT_PREDICTION] = get_top_n_prots(subsequence_aggs)
    json[EXPERIMENT_ENTRY][subsequence_name] = subsequence_dict

#####################################################
#               END "PRIVATE" FUNCTIONS
#####################################################

'''analyze

DESC:
    entry point to save all raw data into an experiment json
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
    files: list of str output files from the scoring algorithm
OPTIONAL:
    predicting_agg_func: str name of the aggregation function to use. Default=sum
    saving_dir: str the name of the directory to save the experiment in. Default=./
RETURNS:
    str file path to the experiment json generated
'''
def analyze(proteins, peptides, files, predicting_agg_func='sum', saving_dir='./'):
    global experiment_json_file_name, experiment_json
    saving_dir = utils.__make_valid_dir_string(saving_dir)
    utils.__make_dir(saving_dir)

    # add header information to the json
    __add_header_info(proteins, peptides, experiment_json)

    # isolate the protein names
    protein_names = [x['name'] for x in proteins]

    for pep in peptides:
        peptide_name = pep['peptide_name']
        pep_related = utils.__get_related_files(files, peptide_name)
        __save_subsequence_info(pep_related, peptide_name, protein_names, experiment_json, predicting_agg_func=predicting_agg_func)
    
    # save to file
    write_raw_json(saving_dir + experiment_json_file_name, experiment_json)
    return(saving_dir + experiment_json_file_name)
