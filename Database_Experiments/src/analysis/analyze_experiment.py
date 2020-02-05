from copy import deepcopy
import utils
from sequence_generation.digest import load_digest
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
    args: dictionary parameters used when running the experiment
    json: dictionary object in which to save the header info
'''
def __add_header_info(proteins, peptides, args, json):
    json[EXPERIMENT_HEADER][EXPERIMENT_PROTEIN_HEADER] = deepcopy(proteins)
    json[EXPERIMENT_HEADER][EXPERIMENT_PEPTIDE_HEADER] = deepcopy(peptides)
    json[EXPERIMENT_HEADER][EXPERIMENT_ARGUMENT_HEADER] = deepcopy(args)

'''__add_subsequence_agg

DESC:
    gets and retrieves all information of the subsequence 
    Updates the global experiment_json object
PARAMS:
    subsequence_files: list of strings filepaths to files of the experiment
    subsequence_name: string name of the subsequence
    protein_names: names of all the proteins 
    json: dictionary object where all items are save
OPTIONAL:
    predicting_agg_func: str name of the function used for aggregation. Default=sum
    mix_in_hybrids: bool whether or not include hybrid proteins in analysis. Default=False
RETURNS: 
    None
'''
def __add_subsequnce_agg(subsequence_files, subsequence_name, protein_names, json, predicting_agg_func='sum', mix_in_hybrids=False):
    agg_func = __z_score_sum if 'z_score_sum' in predicting_agg_func.lower() else ( __product if 'product' in predicting_agg_func.lower() else __sum)
    subsequence_dict = {}
    subsequence_aggs = {}
    for prot_name in protein_names:
        if HYBRID_SEACH_STRING.lower() in str(prot_name).lower() and not mix_in_hybrids:
            continue

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

        ranking[prot] = {}
        for k, scores in kmers.items():
            if k not in tagged_scores:
                tagged_scores[k] = []
            ranking[prot][k] = None
            tagged = [(prot, score, i) for i, score in enumerate(scores)]
            tagged_scores[k] += tagged

    for k, tg in tagged_scores.items():
        tg.sort(reverse=True, key=lambda x: x[1])
        for r, score in enumerate(tg):
            # TODO: make all k-mers with the same score the same rank
            if int(score[2]) == correct_position:
                ranking[score[0]][k] = r + 1 # 1 based score instead of 0
    
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

'''__analyze_subsequence

DESC:
    gets and retrieves all information of the subsequence 
PARAMS:
    pep: dictionary for all information pertaining to the peptide
    files: list of str filepaths for all output files
    protein_names: names of all the proteins 
    json: dictionary object where all items are save
OPTIONAL:
    predicting_agg_func: str name of the function used for aggregation. Default=sum
    mix_in_hybrids: bool whether or not include hybrid proteins in analysis. Default=False
RETURNS: 
    None
'''
def __analyze_subsequence(pep, files, protein_names, json, predicting_agg_func='sum', mix_in_hybrids=False):
    peptide_name = pep['peptide_name']
    parent_name = pep['parent_name']
    if not mix_in_hybrids and HYBRID_SEACH_STRING in parent_name:
        return 
    pep_related = utils.__get_related_files(files, peptide_name)
    __add_subsequnce_agg(pep_related, peptide_name, protein_names, json, predicting_agg_func=predicting_agg_func, mix_in_hybrids=mix_in_hybrids)
    __rank_pep(json, pep)

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
    args: dictionary parameters used when running the experiment
OPTIONAL:
    predicting_agg_func: str name of the aggregation function to use. Default=sum
    saving_dir: str the name of the directory to save the experiment in. Default=./
    mix_in_hybrids: bool whether or not to include using hybrid proteins in analysis. Default=False
    show_all: bool whether or not to show all generated plots. Default=False
    compress: bool to compress subsequence generated directories. Default=True
    hide_hybrids: bool to hide hybrid proteins when plotting. Default=True
RETURNS:
    str file path to the experiment json generated
'''
def analyze(proteins, peptides, files, args, predicting_agg_func='sum', saving_dir='./', mix_in_hybrids=False, show_all=False, compress=True, hide_hybrids=True):
    '''
    1. Perform any aggregations
    2. Rank peptides and do stats
    3. Save to file
    '''
    global experiment_json_file_name, experiment_json
    saving_dir = utils.__make_valid_dir_string(saving_dir)
    utils.__make_dir(saving_dir)

    # add header information to the json
    __add_header_info(proteins, peptides, args, experiment_json)

    # isolate the protein names
    protein_names = [x['name'] for x in proteins]

    # perform analysis on all peptides
    p_counter = 0
    for pep in peptides:
        p_counter += 1
        print('Analyzing peptide {}/{}[{}%]\r'.format(p_counter, len(peptides), int( (float(p_counter)/float(len(peptides))) *100 )), end='')
        __analyze_subsequence(pep, files, protein_names, experiment_json, predicting_agg_func=predicting_agg_func, mix_in_hybrids=mix_in_hybrids)
    
    # save to file
    write_raw_json(saving_dir + experiment_json_file_name, experiment_json)

    # load the experiment and plot it
    plot_experiment(experiment_json, agg_func=predicting_agg_func, show_all=show_all, saving_dir=saving_dir, compress=compress, hide_hybrids=hide_hybrids)

    return(saving_dir + experiment_json_file_name)
