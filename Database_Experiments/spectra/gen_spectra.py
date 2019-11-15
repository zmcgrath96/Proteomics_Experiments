from pyteomics import mass

'''
[{sequence, [], precursor mass}]
'''
def gen_spectra(sequenences):
    spectra = []
    for sequence in sequenences:
        rev_sequence = sequence[::-1]
        this_spectra = []
        this_entry = {'sequence': sequence, 'spectrum': None, 'precursor_mass': None}
        for i in range(1, len(sequence) + 1):
            subsequence = sequence[0:i]
            rev_subsequence = rev_sequence[0:i]
            # charge 1 has an extra proton
            # without a charge and just leaving the ion for b gives us residues
            this_spectra.append(mass.calculate_mass(sequence=subsequence, ion_type='b', charge=1))
            this_spectra.append(mass.calculate_mass(sequence=rev_subsequence, ion_type='y', charge=1))
            this_spectra.append(mass.calculate_mass(sequence=subsequence, ion_type='b', charge=2))
            this_spectra.append(mass.calculate_mass(sequence=rev_subsequence, ion_type='y', charge=2))
        this_spectra.sort()
        this_entry['spectrum'] = this_spectra
        this_entry['precursor_mass'] = mass.calculate_mass(sequence=sequence, charge=2)
        spectra.append(this_entry)

    return spectra