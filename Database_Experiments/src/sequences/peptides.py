from sequences import digest

def gen_peptides(sequence_dict, number_peptides, peptide_index='peptides', min_length=3, save_dir='./'):
    proteins = sequence_dict['sample']['proteins']
    peptides = digest.tryptic(proteins, number_peptides, min_length=min_length, save_dir=save_dir)
    sequence_dict[peptide_index] = peptides
