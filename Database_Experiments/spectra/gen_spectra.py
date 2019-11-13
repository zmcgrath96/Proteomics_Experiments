from pyteomics import mass

def gen_spectra(sequenences):
    spectra = []
    for sequence in sequenences:
        rev_sequence = sequence[::-1]
        this_spectra = []
        for i in range(1, len(sequence) + 1):
            subsequence = sequence[0:i]
            rev_subsequence = rev_sequence[0:i]
            this_spectra.append(round(mass.calculate_mass(sequence=subsequence, ion_type='b'), 5))
            this_spectra.append(round(mass.calculate_mass(sequence=rev_subsequence, ion_type='y'), 5))
        this_spectra.sort()
        spectra.append(this_spectra)

    return spectra