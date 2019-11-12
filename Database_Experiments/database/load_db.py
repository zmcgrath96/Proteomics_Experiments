from pyopenms import FASTAFile

def load_fasta(db_file_path):
    fasta = FASTAFile()
    entries = []
    fasta.load(db_file_path, entries)
    return entries