from sequences import digest

def gen_peptides(sequence_dict, number_peptides, peptide_index='peptides'):
    proteins = sequence_dict['sample']['proteins']
    peptides = digest.tryptic(proteins, number_peptides)
    sequence_dict[peptide_index] = peptides
