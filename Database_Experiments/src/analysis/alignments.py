from copy import deepcopy
from analysis.analysis_utils import get_top_n_prots

###################################################################
#                       CONSTANTS
###################################################################

START_POSITION = 'starting_position'
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
            'starting_position': int,
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

def predict_pep_seqs(prot_scores: list, peptide_dict: dict) -> list:
    '''
    For all proteins try and predict the sequence

    Inputs:
        prot_scores:    list of scores with protein information
        peptide_dict:   dictionary containing protein names (keys) and score information (dict) with kmers (k=keys) and scores (values)
    Outputs:
        list of dictionaries with prediction information
        {
            'starting_position':    int,
            'predicted_length':     int,
            'protein_name':         str,
        }
    '''
    predictions = []
    for tp in prot_scores:
        prot_name = tp[PROTEIN_NAME]
        pos = tp[POSITION]
        prediction = __predict_sequence(peptide_dict[prot_name], pos)
        prediction[PROTEIN_NAME] = prot_name
        predictions.append(prediction)
    return predictions


###################################################################
#                     END PRIVATE FUNCTIONS
###################################################################

def make_sequence_predictions(peptide_dict: dict, agg_func:str, n=5) -> list:
    '''
    make a prediction as to who the parent is and what the sequence is
    
    Inputs:
        peptide_dict:           dictionary containing protein names (keys) and score information (dict) with kmers (k=keys) and scores (values)
        agg_func:               aggregation fucntion name. Used to find the aggregated scores
    kwargs:
        n:                      int top n preditions to make
    Outputs:
        peptide_predition:      list of top n predictions
    '''
    agged = {}
    for p in peptide_dict:
        if p == SAMPLE_PROTEIN_ANALYSIS:
            continue
        agged[p] = deepcopy(peptide_dict[p][agg_func])
    
    # get the best proteins
    top_prots = get_top_n_prots(agged, n=n)

    # now that we have the top proteins, use k-information to try and predict the sequece
    predictions = predict_pep_seqs(top_prots, peptide_dict)

    return predictions

def make_sequence_predictions_ions(peptide_dict: dict, agg_func: str, n=5) -> list:
    '''
    Make best n predictions for both b ions and y ions

    Inputs:
        peptide_dict:           dictionary containing protein names (keys) and score information (dict) with kmers (k=keys) with ions ('b, 'y') and scores (values)
        agg_func:               aggregation fucntion name. Used to find the aggregated scores
    kwargs:
        n:                      int top n preditions to make
    Outputs:
        b_prediction, y_prediction
            both of these are list of top n predictions
    '''
    # separate b and y info
    b_pep_dict = {}
    y_pep_dict = {}
    b_prots = {}
    y_prots = {}
    # {k: {b: , y: }}
    get_ion_from_kmers = lambda x, ion: {k: x[k][ion] for k in x}
    for p in peptide_dict:
        if p == SAMPLE_PROTEIN_ANALYSIS:
            continue
        b_prots[p] = deepcopy(peptide_dict[p][agg_func]['b'])
        y_prots[p] = deepcopy(peptide_dict[p][agg_func]['y'])
        b_pep_dict[p] = get_ion_from_kmers(peptide_dict[p], 'b') 
        y_pep_dict[p] = get_ion_from_kmers(peptide_dict[p], 'y')
    
    # get the best proteins
    top_prots_b = get_top_n_prots(b_prots, n=n)
    top_prots_y = get_top_n_prots(y_prots, n=n)

    # now that we have the top proteins, use k-information to try and predict the sequece
    b_predictions = predict_pep_seqs(top_prots_b, b_pep_dict)
    y_predictions = predict_pep_seqs(top_prots_y, y_pep_dict)

    return b_predictions, y_predictions