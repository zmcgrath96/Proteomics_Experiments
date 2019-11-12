from database import load_db

def load_dbs(iter_of_fasta_files):
    dbs = []
    for fasta_file in iter_of_fasta_files:
        dbs.append(load_db.load_fasta(fasta_file))

    return dbs