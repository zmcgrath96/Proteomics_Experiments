from data.score_utils import __pad_scores
import numpy as np

'''__sum

DESC:
    aggregation function that take the sum of all scores
PARAMS:
    all_scores: iterable of lists of scores (floats)
RETURNS:
    list of floats
'''
def __sum(all_scores):
    score_sum = []
    if type(all_scores) is dict:
        for _, score in all_scores.items():
            score_sum, score = __pad_scores(score_sum, score)
            score_sum = np.add(score_sum, score)
    else:
        for score in all_scores:
            score_sum, score = __pad_scores(score_sum, score)
            score_sum = np.add(score_sum, score)

    return list(score_sum)

'''__product

DESC:
    aggregation function that takes the product of all positive scores
PARAMS:
    all_scores: iterable of lists of scores (floats)
RETURNS:
    list of floats
'''
def __product(all_scores):
    score_prod = []
    if type(all_scores) is dict:
        for _, score in all_scores.items():
            score_prod, score = __pad_scores(score_prod, score, padding=1)
            for i in range(len(score_prod)):
                if score[i] > 0:
                    score_prod *= score
    else:
        for score in all_scores:
            score_prod, score = __pad_scores(score_prod, score, padding=1)
            for i in range(len(score_prod)):
                if score[i] > 0:
                    score_prod *= score

    return list(score_prod)

'''__z_score_sum

DESC:
    aggregation function that takes the z-score-sum of all k-mers
PARAMS:
    all_scores: iterable of lists of scores (floats)
RETURNS:
    list of floats 
'''
def __z_score_sum(all_scores):
    scores = all_scores if not type(all_scores) is dict else [l for _, l in all_scores.items()]
    score_z_sum = []
    s = np.std([item for sublist in scores for item in sublist])
    u = np.mean([item for sublist in scores for item in sublist])
    for scr in scores:
        score_z_sum, scr = __pad_scores(score_z_sum, scr)
        t = [(x-u)/s for x in scr]
        score_z_sum = np.add(score_z_sum, t)
    
    return list(score_z_sum)