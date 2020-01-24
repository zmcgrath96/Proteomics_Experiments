from spectra.gen_spectra import __calc_masses

'''cmp_spectra_spectra

DESC: 
    score two spectra against eachother
    simple additive scoring algorithm
PARAMS:
    spec1: list of floats mass spectra of first sequence
    spec2: list of floats mass spectra of second sequence
RETURNS:
    float score from comparison
'''
def cmp_spectra_spectra(spec1, spec2):
    streak = 0
    last = False
    score = 0
    max_streak = 0
    for mass in spec1:
        if last == True:
            streak += 1
            max_streak = max([streak, max_streak])

        if mass in spec2:
            score += 1
            last = True

        else:
            streak = 0
            last = False
    
    score += max_streak
    divider = min([len(spec1), len(spec2)]) / 2
    return score / divider

'''cmp_string_spectra

DESC:
    compare a string and a spectra together
    uses simple additive scoring
    uses both single and doubly charged ions for mass calculations
PARAMS:
    seq: string sequence of amino acids to convert to mass spectra
    spec: list of floats mass spectra
RETURNS:
    float score from comparison
'''
def cmp_string_spectra(seq, spec):
    spec1 = []
    m1, _ = __calc_masses(seq, 1)
    m2, _ = __calc_masses(seq, 2)
    spec1 = m1 + m2
    return cmp_spectra_spectra(spec1, spec)

'''cmp_string_string

DESC:
    compare the two spectras from two strings
    uses simple additive scoring
PARAMS:
    seq1: string sequence of amino acids
    seq2: string sequence of amino acids
RETURNS:
    float score from comparison
'''
def cmp_string_string(seq1, seq2):
    spec1, spec2 = [], []
    m11, _ = __calc_masses(seq1, 1)
    m12, _ = __calc_masses(seq1, 2)
    m21, _ = __calc_masses(seq2, 1)
    m22, _ = __calc_masses(seq2, 2)
    spec1 = m11 + m12 
    spec2 = m21 + m22 
    return cmp_spectra_spectra(spec1, spec2)

