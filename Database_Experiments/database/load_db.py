from pyteomics import fasta

def load_fasta(db_file_path):
    entries = []
    with fasta.read(db_file_path) as db:
        for entry in db:
            entries.append(entry)
    return entries