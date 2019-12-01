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
    func = average
    if 'sum' in measure.lower():
        func = sum 
    elif 'max' in measure.lower():
        func = max
    
    l = [0]
    v = None
    if type(scores) is dict:
        for key in scores:
            if func(scores[key]) > func(l):
                l = scores[key]
                v = key

    elif type(scores) is list:
        for i, score in enumerate(scores):
            if func(score) > func(l):
                l = score 
                v = i
    else:
        l = []
        v = 0

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
        