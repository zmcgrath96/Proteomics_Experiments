from file_io import mzML, fasta, csv
from scoring.compare import cmp_string_spectra


def search_proteins(spectrum: dict, database: list) -> dict:
    '''
    Find the highest scoring protein from a spectrum in a database

    Inputs:
        spectrum:   dictionary of the form {'level': , 'scan_no': , 'spectrum': }
        database:   list of dictionaries read from a fasta file of form
                    {
                        'identifier': any,
                        'sequence': str, 
                        'name': str,
                        'scan_no': int,
                        'level': int
                    }
    Outputs:
        Dictionary of the hightest match of the form
        {
            'protein_id': any,
            'protein_name': str,
            'protein_sequence': str,
            'ms_level': int, 
            'scan_no': int,
            'score': float
        }
    '''
    best_match = {}
    for prot in database:
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

def search_database(spec_file: str, database: list, output_name: str) -> str:
    '''
    Find the hightest scoring subsequence for each spectra
    
    Inputs:
        spec_file:      str path to .mzML file
        database:       list of dictionaries read from a fasta file of form
                    {
                        'identifier': any,
                        'sequence': str, 
                        'name': str,
                        'scan_no': int,
                        'level': int
                    }
        output_name:    string name of file to save with path 
    Outputs: 
        string path to the saved file containing scores
    '''
    specs = mzML.read(spec_file)
    results = []
    for spec in specs:
        results.append(search_proteins(spec, database))

    return csv.write_iter_of_dicts(results, output_name)
