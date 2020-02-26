import pyopenms
import os
import shutil
from subprocess import call
from utils import __make_dir, __make_valid_dir_string, __is_gzipped, __gunzip, __gzip
from scoring import search

crux_to_rm = ['tide-search.decoy.txt', 'tide-search.log.txt', 'tide-search.params.txt']
scoring_functions = ['crux', 'custom']

#######################################################################################
#                       BEGIN "PRIVATE" FUNCTIONS
#######################################################################################

def __parse_spectrum_name(spec_name):
    return str(spec_name.split('/')[-1]).lower().replace('.mzml', '')

def __parse_db_name(db_name):
    return  str(db_name.split('/')[-1]).replace('.fasta', '')

def __index_db_files(path_to_crux_cmd, db_files):
    idx_names = []
    num_dbs = len(db_files)

    for i, db_file in enumerate(db_files):
        print('On database {}/{}[{}%]\r'.format(i+1, num_dbs, int(((i+1)/num_dbs) * 100)), end="")
        
        this_output_dir = '/'.join(str(db_file).split('/')[:-1])
        this_output_dir = __make_valid_dir_string(this_output_dir) + 'indexed/'
        __make_dir(this_output_dir)
        idx_name = this_output_dir + str(db_file).replace('.fasta', '_index').split('/')[-1]

        indx_cmd = [
            path_to_crux_cmd, 
            'tide-index', 
            db_file, 
            idx_name, 
            '--min-length', '2', 
            '--min-mass', '50', 
            '--output-dir', this_output_dir, 
            '--overwrite', 'T', 
            '--min-peaks', '2', 
            '--precursor-window', '1000000000',
            '--enzyme', 'no-enzyme', 
            '--verbosity', '0'
        ]
        call(indx_cmd)
        # remove extra output files
        os.remove(this_output_dir + 'tide-index.params.txt')
        os.remove(this_output_dir + 'tide-index.log.txt')

        idx_names.append(idx_name)
    return idx_names

def __remove_indices(index_file):
    if isinstance(index_file, list):
        index_file = index_file[0]
    rm_dir = '/'.join(index_file.split('/')[:-1])
    shutil.rmtree(rm_dir)

'''crux_search

DESC:
    use the crux tool to score spectra against databases
Inputs:
    spectra_files: list of str paths to all the spectra (.mzML) files
    database_files: list of str paths to all the database (.fasta) files
    path_to_crux_cmd: str path to the executable for crux.
    output_dir: str path to the directory to save files
kwargs:
    compress: bool compress the output result. Default=True
Outputs:
    list of str of output files
'''
def __crux_search(spectra_files, database_files, path_to_crux_cmd, output_dir, compress=True):
    output_dir = __make_valid_dir_string(output_dir) + 'search_output/'
    __make_dir(output_dir)
    spec_dir = '/'.join(spectra_files[0].split('/')[:-1])

    is_compressed = __is_gzipped(spectra_files[0])
    print('Pre-indexing database files...')
    indexed_db_files = __index_db_files(path_to_crux_cmd, database_files)
    print('\nDone. Scoring..')

    output_count = 0
    output_files = []
    num_dbs = len(indexed_db_files)
    num_specs = len(spectra_files)

    for i, spec_file in enumerate(spectra_files):
        spec_file = spec_file if not is_compressed else __gunzip(spec_file)
        for j, database_file in enumerate(indexed_db_files):
            this_db_name = __parse_db_name(database_file)
            print('On spectrum {}/{}[{}%]\tOn database file {}/{}[{}%]\r'.format(i+1, num_specs, int(((i+1)/num_specs) * 100), j+1, num_dbs, int(((j+1)/num_dbs)*100)), end="")
            this_output_dir = output_dir + '{}_vs_{}'.format(__parse_spectrum_name(spec_file), this_db_name)
            search_cmd = [
                path_to_crux_cmd, 
                'tide-search', 
                spec_file, 
                database_file, 
                '--min-length', '2', 
                '--min-mass', '50', 
                '--output-dir', this_output_dir, 
                '--overwrite', 'T', 
                '--min-peaks', '2', 
                '--precursor-window', '1000000000',
                '--enzyme', 'no-enzyme', 
                '--verbosity', '0'
                ] 
            call(search_cmd)
            output_count += 1
            o = this_output_dir + '/tide-search.target.txt' if not compress else __gzip(this_output_dir + '/tide-search.target.txt')
            o_tsv = o.replace('.txt', '.tsv')
            os.rename(o, o_tsv)
            output_files.append(o_tsv)

            # saving space, so should remove the extra stuff
            if is_compressed:
                for rm in crux_to_rm:
                    os.remove(this_output_dir + '/' + rm)

        # if the files were compressed, we were trying to save disk space so just remove the mzml files
        is_compressed and os.remove(spec_file)

    is_compressed and os.rmdir(spec_dir)
    __remove_indices(indexed_db_files)
    return output_files

'''__custom_search

DESC:
    custom search on peptides
Inputs:
    spectra_files: list of strings paths to all .mzML files
    database_files: list of strings paths to all .fasta files
    output_dir: str path to directory to save results under
OPTONAL:
    compress: bool whether or not to compress results. Default=True
Outputs:
    list of str of output files
'''
def __custom_search(spectra_files, database_files, output_dir, compress=True):
    output_dir = __make_valid_dir_string(output_dir)
    output_files = []

    no_spec = len(spectra_files)
    no_db = len(database_files)

    for spec_no, spectra_file in enumerate(spectra_files):
        spectra_file = spectra_file if not __is_gzipped(spectra_file) else __gunzip(spectra_file)

        for db_no, database_file in enumerate(database_files):
            print('On spectrum: {}/{} [{}%]   On database: {}/{} [{}%]\r'.format(spec_no, no_spec, int(float(spec_no) / float(no_spec) * 100), db_no, no_db, int(float(db_no) / float(no_db) * 100)), end='')

            output_name = output_dir + 'search_output/' + '{}_vs_{}'.format(__parse_spectrum_name(spectra_file), __parse_db_name(database_file))
            database_file = database_file if not __is_gzipped(database_file) else __gunzip(database_file)
            output_file = search.search_files(spectra_file, database_file, output_name)
            output_file = output_file if not compress else __gzip(output_file)
            output_files.append(output_file)
    return output_files
#######################################################################################
#                        END "PRIVATE" FUNCTIONS
#######################################################################################

'''score_peptides

DESC:
    use the crux tool to score spectra against databases
Inputs:
    spectra_files: list of str paths to all the spectra (.mzML) files
    database_files: list of str paths to all the database (.fasta) files
    output_dir: str path to the directory to save files
kwargs:
    compress: bool compress the output result. Default=True
    score_func: str determine which scoring function to use. Default='custom'
    path_to_crux_cmd: str path to the executable for crux. Default=''
Outputs:
    list of str of output files
'''
def score_peptides(spectra_files, database_files, output_dir, compress=True, score_func='custom', path_to_crux_cmd=''):
    score_func = 'custom' if score_func.lower() not in scoring_functions else score_func.lower()
    if score_func == 'crux':
        return __crux_search(spectra_files, database_files, path_to_crux_cmd, output_dir, compress=compress)
    else:
        return __custom_search(spectra_files, database_files, output_dir, compress=compress)