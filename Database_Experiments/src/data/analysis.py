from numpy import argmax, average
from math import inf

'''__get_argmax_max

DESC:
    find the argmax and the max of an iterable
PARAMS:
    i: iterable to return max and argmax on
RETURNS:
    int, val: argmax, max
'''
def __get_argmax_max(i):
    return argmax(i), max(i)

'''get_highest_scoring

DESC:
    return the list and key of the best scoring list of scores
PARAMS:
    scores: iterable of lists each entry of the lists should be a number. If scores is not an iterable, [], 0 is returned
OPTIONAL:
    measure: string determines which metric is used for the calculate best score. Can be 'max', 'average', 'sum'. Default=average
RETURNS:
    list, val: the list and the index or key of the best score list
'''
def get_highest_scoring(scores, measure='average'):
    func = sum if 'sum' in measure.lower() else (max if 'max' in measure.lower() else average)
    
    l = None
    v = None
    if type(scores) is dict:
        for key, score in scores.items():
            if l is None:
                l = score
                v = key
            if func(score) > func(l):
                l = score
                v = key
    else:
        for i, score in enumerate(scores):
            if l is None:
                l = score
                v = i
            if func(score) > func(l):
                l = score 
                v = i
    return l, v

'''get_top_n

DESC:
    returns the top n list and keys or indices of an interable of lists
PARAMS:
    scores: iterable of lists each entry of the lists should be a number. If scores is not an iterable, [] is returned
OPTIONAL:
    measure: string determines which metric is used for the calculate best score. Can be 'max', 'average', 'sum'. Default=average
    n: int number of top scores to return. Default=5
RETURNS:
    [(list, val)]: list of tuples with the list, index or key of the top n in order
'''
def get_top_n(scores, measure='average', n=5):
    top_n = []
    n = len(scores) if len(scores) < n else n

    for _ in range(n):
        l, v = get_highest_scoring(scores, measure=measure)
        t = (l, v)
        top_n.append(t)
        del scores[v]
    
    return top_n
        