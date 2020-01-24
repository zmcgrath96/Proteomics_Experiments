from file_io import mzML, fasta, csv
from scoring.score import cmp_string_spectra

'''search_proteins

DESC:
    find the highest scoring protein from a spectrum in a database
PARAMS:
    spectrum: dictionary of the form {'level': , 'scan_no': , 'spectrum': }
    database: .fasta file of proteins to search through
RETURNS:

'''
def search_proteins(spectrum, database):
    db = fasta.read(database)
    best_match = {}
    for prot in db:
        score = cmp_string_spectra(prot['sequence'], spectrum['spectrum'])
        if best_match == {} or best_match['score'] < score:
            best_match = {
                'protein_id': prot['identifier'],
                'protein_name': prot['name'], 
                'protein_sequence': prot['sequence'], 
                'ms_level': spectrum['level'],
                'scan_no': spectrum['scan_no'], 
                'score': score
            }

    return best_match


'''search_files

DESC:
    find the hightest scoring subsequence for each spectra
PARAMS:
    spec_file: str path to .mzML file
    database_file: str path to .fasta file
OPTONAL:
    output_dir: str path to directory to save output files to. Default='./'
    compress: bool to compress output file or not. Default=True
RETURNS: 
    str path to output file
'''
def search_files(spec_file, database_file, output_name):
    specs = mzML.read(spec_file)
    results = []
    for spec in specs:
        results.append(search_proteins(spec, database_file))

    return csv.write_iter_of_dicts(results, output_name)
