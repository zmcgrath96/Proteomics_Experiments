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

def __b_ions(sequence: str, charge=None): 
    '''
    Calculate the masses for b_ions
    
    Inputs:
        sequence:  string amino acid sequence to calculate
    kwargs:
        charge:    int charge to calculate. Possible types are {1, 2}. Default is both
    Outputs:
        list of floats
    '''
    masses = []
    length = len(sequence)
    
    if charge is None or charge == 1:
        #b+
        total = 1.007825035 - 0.0005486 #for the H to turn the residue NH on the N-terminus into NH2
        for i in range (0, length):
            total += amino_acids[sequence[i]]
            masses.append(total)
            #Since z (the charge) is equal to one, the total here is the m/z
            
    if charge is None or charge == 2:
        #b++
        total = 2 * 1.007825035 - 2 * 0.0005486 #adding one more proton this time to make it doubly charged
        for i in range (0, length):
            total += amino_acids[sequence[i]]
            masses.append(total/2)
            
    return masses

def __y_ions(sequence: str, charge=None): 
    '''
    Calculate the masses for y_ions
    
    Inputs:
        sequence:  string amino acid sequence to calculate
    kwargs:
        charge:    int charge to calculate. Possible types are {1, 2}. Default is both
    Outputs:
        list of floats
    '''
    masses = []
    length = len(sequence)
    
    if charge is None or charge == 1:
        #y+
        total = 3 * 1.007825035 + 15.99491463 - 0.0005486 #for the OH to turn the residue CO on the C-terminus into COOH + 1 proton to make NH into NH2 and 1 proton make positively charged
        for i in range (0,length):
            total += amino_acids[sequence[length-i-1]]
            masses.append(total)
            
    if charge is None or charge == 2:
        #y++
        total = 4 * 1.007825035 + 15.99491463 - 2 * 0.0005486 #another proton to make doubly charged
        for i in range (0, length):
            total += amino_acids[sequence[length-i-1]]
            masses.append(total/2)
            
    return masses

def calc_masses(sequence: str, charge=None, ion=None) -> (list, float):
    '''
    Calculate the molecular weight (Da) of an Amino Acid sequence
    
    Inputs:
        sequence:   string amino acid sequence to calculate mass of
    kwargs:
        charge:     int charge to calculate. Possible values are {1, 2}. Default is both
        ion:        string ion type to calculate. Possible values are {'b', 'y'}. Default is both
    Output:
        (masses, precursor_mass)
        masses:         list of floats of spectrum calculated
        precursor_mass: float precursor mass of the entire amino acid sequence
    '''
    masses = []

    length = len(sequence)
    total = 2 * 1.007825035 + 15.99491463 #This is the mass of water. Adding the mass of water to the sum of all the residue masses gives the mass of the peptide.
    for i in range(length):
        total +=  amino_acids[sequence[i]]

    pre_mz_charge = 2 if charge is None else charge
    pre_mz = (total+pre_mz_charge*1.0072764)/pre_mz_charge   
    
    if ion is None or ion == 'b': 
        masses += __b_ions(sequence, charge=charge)
        
    if ion is None or ion == 'y': 
        masses += __y_ions(sequence, charge=charge)
        
    return masses, pre_mz

def gen_spectrum(sequence: str, charge=None, ion=None) -> list:
    '''
    Generate a spectrum for a single sequence. Includes singly and doubly charged masses
    
    Inputs:
        sequence: string amino acid sequence to calculate spectra for
    kwargs:
        charge:   int charge value to calculate masses for. Possible types are {1, 2}. Default is both
        ion:      string ion type to calculate masses for. Possible types are {'b', 'y'}. Default is both
    Outputs:
        dictionary with the following values 
        {
            'spectrum': list of floats,
            'precursor_mass': float,
        }
    '''
    
    this_entry = {}
    masses, pre_mz = calc_masses(sequence, charge=charge, ion=ion)
    masses.sort()
    this_entry['spectrum'] = masses
    this_entry['precursor_mass'] = pre_mz
    return this_entry

def gen_spectra(sequences: list, charge=None, ion=None) -> list:
    '''
    Generates mass spectra for a list of sequences. Includes singly and doubly charged masses

    Inputs:
        sequences: list of strings sequences to generate spectra for
    kwargs:
        charge:   int charge value to calculate masses for. Possible types are {1, 2}. Default is both
        ion:      string ion type to calculate masses for. Possible types are {'b', 'y'}. Default is both
    Outputs:
        list of dictionaries of the form {'spectrum': list of floats, 'precursor_mass': float}
    '''
    return [gen_spectrum(seq, charge=charge, ion=ion) for seq in sequences]