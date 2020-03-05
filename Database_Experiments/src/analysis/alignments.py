from copy import deepcopy
from analysis.analysis_utils import get_top_n_prots

###################################################################
#                       CONSTANTS
###################################################################

START_POSITION = 'starting_pos'
PREDICTED_LENGTH = 'predicted_length'
SAMPLE_PROTEIN_ANALYSIS = 'analysis'
PROTEIN_NAME = 'protein_name'
POSITION = 'position'

###################################################################
#                     END CONSTANTS
###################################################################

###################################################################
#                       PRIVATE FUNCTIONS
###################################################################

def __predict_sequence(prot_info: dict, starting_pos: int) -> dict:
    '''
    make a prediction on what the peptide sequence is

    Inputs:
        prot_info:      dictionary with the kmer scoring info
        starting_pos:   int position of the higest score
    Outputs:
        dictionary with the prediction. 
        {
            'starting_pos': int,
            'predicted_length': int
        }
    '''
    max_score_k = 0
    max_score = -10
    get_k = lambda sk: int(sk[sk.index('=')+1:])
    for ke in prot_info:
        if '=' not in ke:
            continue
        # check to see that the starting position is in the length. If its not, we know the peptide is shorter than that k
        if starting_pos >= len(prot_info[ke]):
            continue
        k = get_k(ke)
        # keep the current champion score if the max score is at least equal to the current score
        max_score, max_score_k = (max_score, max_score_k) if max_score >= prot_info[ke][starting_pos] else (prot_info[ke][starting_pos], k)
    # we now have the k where the score peaked, so we should be able to make dumb prediction off this
    return {START_POSITION: starting_pos, PREDICTED_LENGTH: max_score_k}

###################################################################
#                     END PRIVATE FUNCTIONS
###################################################################

def make_sequence_predictions(peptide_dict: dict, agg_fucn:str, n=5) -> dict:
    '''
    make a prediction as to who the parent is and what the sequence is
    
    Inputs:
        peptide_dict:           dictionary containing protein names (keys) and score information (dict) with kmers (k=keys) and scores (values)
        agg_func:               aggregation fucntion name. Used to find the aggregated scores
    kwargs:
        n:                      int top n preditions to make
    Outputs:
        peptide_predition:      dictionary of top n predictions
    '''
    agged = {}
    for p in peptide_dict:
        if p == SAMPLE_PROTEIN_ANALYSIS:
            continue
        agged[p] = deepcopy(peptide_dict[p][agg_fucn])
    
    # get the best proteins
    top_prots = get_top_n_prots(agged, n=n)

    # now that we have the top proteins, use k-information to try and predict the sequece
    predictions = []
    for tp in top_prots:
        prot_name = tp[PROTEIN_NAME]
        pos = tp[POSITION]
        prediction = __predict_sequence(peptide_dict[prot_name], pos)
        prediction[PROTEIN_NAME] = prot_name
        predictions.append(prediction)

    return predictions