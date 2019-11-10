import json

def gen_spectra(masses_file, sequenences):
    masses = None
    with open(masses_file, 'r') as mf:
        masses = json.load(mf)

    spectra = []

    for sequenence in sequenences:
        reversed_seq = sequenence[::-1]
        forward_masses = []
        reverse_masses = []
        current_mass = 0

        for aa in sequenence:
            current_mass += masses[aa]
            forward_masses.append(current_mass)

        current_mass = 0
        for aa in reversed_seq:
            current_mass += masses[aa]
            reverse_masses.append(current_mass)

        all_masses = forward_masses + reverse_masses
        all_masses.sort()
        spectra.append(all_masses[:-1])

    return spectra