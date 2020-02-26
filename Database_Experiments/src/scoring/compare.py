from spectra.gen_spectra import calc_masses

'''cmp_spectra_spectra

DESC: 
    score two spectra against eachother
    simple additive scoring algorithm with a length divider
PARAMS:
    spec: list of floats mass spectra of first sequence
    reference: list of floats mass spectra of second sequence
RETURNS:
    float score from comparison
'''
def cmp_spectra_spectra__JAN_2020(spec, reference):
    if len(spec) == 0 or len(reference) == 0:
        return
    streak = 0
    last = False
    score = 0
    max_streak = 0
    for mass in spec:
        if last == True:
            streak += 1
            max_streak = max([streak, max_streak])

        if mass in reference:
            score += 1
            last = True

        else:
            streak = 0
            last = False
    
    score += max_streak
    divider = min([len(spec), len(reference)]) / 2
    return score / divider


def cmp_spectra_spectra(spec: list, reference: list) -> float:
    '''
    CREATED FEB 26 2020
    Score two spectra against eachother. Simple additive scoring with bonuses for streaks

    Inputs:
        spec:       list of floats (from mass spectra)
        reference:  list of floats (calculated from protein sequence)
    Outputs:
        score:      float score 
    '''
    if len(spec) == 0 or len(reference) == 0:
        return
    streak = 0
    last = False
    score = 0
    max_streak = 0
    for mass in spec:
        if last == True:
            streak += 1
            max_streak = max([streak, max_streak])

        if mass in reference:
            score += 1
            last = True

        else:
            streak = 0
            last = False
    
    score += max_streak
    return score 


'''cmp_string_spectra

DESC:
    compare a string and a spectra together
    uses simple additive scoring
    uses both single and doubly charged ions for mass calculations
PARAMS:
    seq: string sequence of amino acids to convert to mass spectra
    ref_spec: list of floats mass spectra
RETURNS:
    float score from comparison
'''
def cmp_string_spectra(seq, ref_spec):
    spec1 = []
    m1, _ = calc_masses(seq, 1)
    m2, _ = calc_masses(seq, 2)
    spec1 = m1 + m2
    return cmp_spectra_spectra(spec1, ref_spec)

'''cmp_string_string

DESC:
    compare the two spectras from two strings
    uses simple additive scoring
PARAMS:
    seq: string sequence of amino acids
    ref_seq: string sequence of amino acids
RETURNS:
    float score from comparison
'''
def cmp_string_string(seq, ref_seq):
    spec1, spec2 = [], []
    m11, _ = calc_masses(seq, 1)
    m12, _ = calc_masses(seq, 2)
    m21, _ = calc_masses(ref_seq, 1)
    m22, _ = calc_masses(ref_seq, 2)
    spec1 = m11 + m12 
    spec2 = m21 + m22 
    return cmp_spectra_spectra(spec1, spec2)

