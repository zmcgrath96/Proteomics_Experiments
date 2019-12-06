from copy import deepcopy
import utils
from sequences.digest import load_digest
from data.write_output import write_raw_json
from data import score_utils

#######################################################
#                   CONSTANTS
#######################################################
experiment_json_file_name = 'experiment_data.json'
EXPERIMENT_ENTRY = 'experiment'
EXPERIMENT_HEADER = 'header'
EXPERIMENT_PROTEIN_HEADER = 'proteins'
EXPERIMENT_PEPTIDE_HEADER = 'peptides'
SAMPLE_ENTRY = 'sample'
SAMPLE_PROTEINS = 'proteins'
SAMPLE_HYBRID_ENTRY = 'hybrid'
SAMPLE_HYBRID_SEQUENCE = 'sequence'
SAMPLE_HYBRID_PARENT = 'full_parent'
SAMPLE_HYBRID_INDICES = 'parent_indices'
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

'''__parse_output_name

DESC:
    get a usable label for a comparision 
PARAMS: 
    crux_output: a string that will use the parent dir of the output as the label
RETURNS:
    a string with the title 
'''
def __parse_output_name(crux_output):
    return crux_output.split('/')[-2]

'''__get_peptide
DESC:
    given a name, get the entry associated with it
PARAMS:
    peptide_name: string name of a peptide entry (peptide_{integer})
RETURNS:
    dictionary entry if it exists else None
'''
def __get_peptide(peptide_name, saving_dir='./'):
    global digests
    saving_dir = utils.__make_valid_dir_string(saving_dir)
    if not digests:
        file_name =  saving_dir + 'digestion.tsv'
        digests = load_digest(file_name)

    return digests[peptide_name] if peptide_name in digests else None

'''__add_header_info
DESC:
    add all header info to the experiment json
'''
def __add_header_info(sequences, saving_dir='./'):
    global digests
    experiment_json[EXPERIMENT_HEADER][EXPERIMENT_PROTEIN_HEADER] = deepcopy(sequences[SAMPLE_ENTRY][SAMPLE_PROTEINS])
    _ = __get_peptide('', saving_dir=saving_dir)
    for key in digests:
        experiment_json[EXPERIMENT_HEADER][EXPERIMENT_PEPTIDE_HEADER].append(deepcopy(digests[key]))

    hybrid_entry = sequences[SAMPLE_HYBRID_ENTRY]
    experiment_json[EXPERIMENT_HEADER][EXPERIMENT_PEPTIDE_HEADER].append({
        "peptide_name": "hybrid",
        "peptide_sequence": hybrid_entry[SAMPLE_HYBRID_SEQUENCE], 
        "parent_name": "hybrid_parent",
        "parent_sequence": hybrid_entry[SAMPLE_HYBRID_PARENT],
        "start_index": hybrid_entry[SAMPLE_HYBRID_INDICES]["start"],
        "end_index": hybrid_entry[SAMPLE_HYBRID_INDICES]["end"]
        })

'''__save_subsequence_info

DESC:
    gets and retrieves all information of the subsequence 
    Updates the global experiment_json object
PARAMS:
    subsequence_files: list of strings filepaths to files of the experiment
    subsequence_name: string name of the subsequence
    protein_names: names of all the proteins 
RETURNS: 
    None
'''
def __save_subsequence_info(subsequence_files, subsequence_name, protein_names):
    global experiment_json
    subsequence_dict = {}
    for prot_name in protein_names:
        subsequence_dict[prot_name] = {}
        prot_with_subseq = utils.__get_related_files(subsequence_files, prot_name)
        for f in prot_with_subseq:
            scores, _, _ = score_utils.__get_scores_scan_pos_label(f)
            mer = str([int(j) for j in str(f).split('_') if j.isdigit()][0])
            k = 'k=' + mer
            subsequence_dict[prot_name][k] = scores
    experiment_json[EXPERIMENT_ENTRY][subsequence_name] = subsequence_dict

#####################################################
#               END "PRIVATE" FUNCTIONS
#####################################################

'''save

DESC:
    entry point to save all raw data into an experiment json
PARAMS:
    experiment, files, protein_names, subsequence_prefix, num_subsequences, hybrid_prefix, sequences, agg_func='sum', show_all=False, saving_dir='./', use_top_n=False, n=5, measure='average'
RETURNS:
    Path to the experiment json
'''
def save(experiment, files, protein_names, subsequence_prefix, num_subsequences, hybrid_prefix, sequences, saving_dir='/'):
    global experiment_json_file_name, experiment_json
    #create the saving directory
    saving_dir = utils.__make_valid_dir_string(saving_dir)
    utils.__make_dir(saving_dir)
    __add_header_info(sequences, saving_dir=saving_dir)

    # do everything for the hybrid sequence first
    hybrid_related = utils.__get_related_files(files, hybrid_prefix)
    __save_subsequence_info(hybrid_related, hybrid_prefix, protein_names)

    # do the rest of the peptides
    peptide_names = []
    for x in range(num_subsequences):
        num = str(x) if x > 9 else '0' + str(x)
        peptide_names.append('{}_{}'.format(subsequence_prefix, num))
    for peptide_name in peptide_names:
        pep_related = utils.__get_related_files(files, peptide_name)
        __save_subsequence_info(pep_related, peptide_name, protein_names)
    
    # save to file
    write_raw_json(saving_dir + experiment_json_file_name, experiment_json)
    return(saving_dir + experiment_json_file_name)
