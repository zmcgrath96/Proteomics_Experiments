from pyteomics import mass

def gen_spectra(sequenences):
    spectra = []
    for sequence in sequenences:
        rev_sequence = sequence[::-1]
        this_spectra = []
        for i in range(1, len(sequence) + 1):
            this_spectra.append(mass.calculate_mass(spectra=sequence[0:i]))
            this_spectra.append(mass.calculate_mass(spectra=rev_sequence[0:i]))
        spectra.append(this_spectra)

    return spectra