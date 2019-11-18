import pyopenms
import os
from subprocess import call

def score_peptides(spectra_files, database_files, path_to_crux_cmd, output_dir):
    # str_spectra_files = ' '.join(spectra_files)
    output_count = 0
    #index_prefix = 'index_'
    for database_file in database_files:
        for spec_file in spectra_files:
            #curr_index = index_prefix + str(output_count)
            #index_cmd = [path_to_crux_cmd, 'tide-index', database_file, curr_index, '--min-length', '3', '--overwrite', 'T', '--min-mass', '100'] #'{} tide-index {} {} --min-length 3 --overwrite T'.format(path_to_crux_cmd, database_file, curr_index)
            search_cmd = [
                path_to_crux_cmd, 
                'tide-search', 
                spec_file, 
                database_file, 
                '--min-length', '3', 
                '--min-mass', '100', 
                '--output-dir', '{}_{}'.format(output_dir, output_count), 
                '--overwrite', 'T', 
                '--min-peaks', '2', 
                '--precursor-window', '1000'
                ] #'{} tide-search {} {} --output-dir {}_{} --overwrite T --min-peaks 2'.format(path_to_crux_cmd, str_spectra_files, database_file, output_dir, output_count)
            #call(index_cmd)
            call(search_cmd)
            output_count += 1