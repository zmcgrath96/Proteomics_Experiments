from pyteomics import mass

amino_acids={
    "A":71.037114,
    "R":156.101111,
    "N":114.042927,
    "D":115.026943,
    "C":103.009185,
    "E":129.042593,
    "Q":128.058578,
    "G":57.021464,
    "H":137.058912,
    "I":113.084064,
    "L":113.084064,
    "K":128.094963,
    "M":131.040485,
    "F":147.068414,
    "P":97.052764,
    "S":87.032028,
    "T":101.047679,
    "U":150.95363,
    "W":186.079313,
    "Y":163.06332,
    "V":99.068414
}

'''__calc_masses

DESC:
    calculates the masses/spectrum for a sequence
PARAMS:
    sequence: str amino acid sequence to change to list of masses
    charge: int charge value to calculate masses for
RETURNS:
    list of floats, float       spectrum and the precursor mass 
'''
def __calc_masses(sequence, charge):
    masses = []

    length = len(sequence)
    total = 2 * 1.007825035 + 15.99491463 #This is the mass of water. Adding the mass of water to the sum of all the residue masses gives the mass of the peptide.
    for i in range(length):
        total +=  amino_acids[sequence[i]]

    pre_mz = (total+charge*1.0072764)/charge   

    if charge == 1:
        #b+
        total = 1.007825035 - 0.0005486 #for the H to turn the residue NH on the N-terminus into NH2
        for i in range (0, length):
            total += amino_acids[sequence[i]]
            masses.append(total)
            #Since z (the charge) is equal to one, the total here is the m/z

        #y+
        total = 3 * 1.007825035 + 15.99491463 - 0.0005486 #for the OH to turn the residue CO on the C-terminus into COOH + 1 proton to make NH into NH2 and 1 proton make positively charged
        for i in range (0,length):
            total += amino_acids[sequence[length-i-1]]
            masses.append(total)

    elif charge == 2:
        #b++
        total = 2 * 1.007825035 - 2 * 0.0005486 #adding one more proton this time to make it doubly charged
        for i in range (0, length):
            total += amino_acids[sequence[i]]
            masses.append(total/2)

        #y++
        total = 4 * 1.007825035 + 15.99491463 - 2 * 0.0005486 #another proton to make doubly charged
        for i in range (0, length):
            total += amino_acids[sequence[length-i-1]]
            masses.append(total/2)
        #The masses you get exactly match Spectrum Mill. To get this, I had to make sure to use the mass of H+ and the mass of H when appropriate.

    return masses, pre_mz


# def __calc_masses_old(sequence):
#     rev_sequence = sequence[::-1]
#     this_spectra = []
#     this_entry = {'sequence': sequence, 'spectrum': None, 'precursor_mass': None}
#     for i in range(1, len(sequence) + 1):
#         subsequence = sequence[0:i]
#         rev_subsequence = rev_sequence[0:i]
#         # charge 1 has an extra proton
#         # without a charge and just leaving the ion for b gives us residues
#         this_spectra.append(mass.calculate_mass(sequence=subsequence, ion_type='b', charge=1))
#         this_spectra.append(mass.calculate_mass(sequence=rev_subsequence, ion_type='y', charge=1))
#         this_spectra.append(mass.calculate_mass(sequence=subsequence, ion_type='b', charge=2))
#         this_spectra.append(mass.calculate_mass(sequence=rev_subsequence, ion_type='y', charge=2))
#     this_spectra.sort()
#     this_entry['spectrum'] = this_spectra
#     this_entry['precursor_mass'] = mass.calculate_mass(sequence=sequence, charge=2)
#     return this_spectra

'''gen_spectra

DESC:
    generates mass spectra for sequences
PARAMS:
    sequences: list of strings sequences to generate spectra for
RETURNS:
    list of dictionaries of the form {'spectrum': list of floats, 'precursor_mass': float}
'''
def gen_spectra(sequenences):
    spectra = []
    for sequence in sequenences:
        this_entry = {}
        mass_1, _ = __calc_masses(sequence, 1)
        mass_2, pre_mz = __calc_masses(sequence, 2)
        this_spectra = mass_1 + mass_2
        this_spectra.sort()
        this_entry['spectrum'] = this_spectra
        this_entry['precursor_mass'] = pre_mz
        spectra.append(this_entry)

    return spectra