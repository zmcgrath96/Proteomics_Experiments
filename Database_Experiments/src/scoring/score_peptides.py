import pyopenms
import os
from subprocess import call
from utils import __make_dir, __make_valid_dir_string, __is_gzipped, __gunzip, __gzip

def __parse_spectrum_name(spec_name):
    return str(spec_name.split('/')[-1]).lower().replace('.mzml', '')

def __parse_db_name(db_name):
    return  str(db_name.split('/')[-1]).replace('.fasta', '')

def score_peptides(spectra_files, database_files, path_to_crux_cmd, output_dir):
    output_dir = __make_valid_dir_string(output_dir) + 'search_output/'
    __make_dir(output_dir)

    is_compressed = __is_gzipped(spectra_files[0])

    output_count = 0
    output_files = []
    num_dbs = len(database_files)
    num_specs = len(spectra_files)
    for i, database_file in enumerate(database_files):
        this_db_name = __parse_db_name(database_file)
        for j, spec_file in enumerate(spectra_files):
            print('On database file {}/{}[{}%]\tOn spectrum {}/{}[{}%]\r'.format(i+1, num_dbs, int(((i+1)/num_dbs) * 100), j+1, num_specs, int(((j+1)/num_specs)*100)), end="")
            this_output_dir = output_dir + '{}_vs_{}'.format(__parse_spectrum_name(spec_file), this_db_name)
            spec_file = spec_file if not is_compressed else __gunzip(spec_file)
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
            output_files.append(this_output_dir + '/tide-search.target.txt')
            is_compressed and __gzip(spec_file)

    return output_files