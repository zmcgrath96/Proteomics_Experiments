import pyopenms
import os

def score_peptides(spectra_files, database_files, path_to_crux_cmd, output_dir):
    str_spectra_files = ' '.join(spectra_files)
    output_count = 0
    index_prefix = 'index_'
    for database_file in database_files:
        curr_index = index_prefix + str(output_count)
        index_cmd = '{} tide-index {} {} --min-length 3 --overwrite T'.format(path_to_crux_cmd, database_file, curr_index)
        search_cmd = '{} tide-search {} {} --output-dir {}_{} --overwrite T'.format(path_to_crux_cmd, str_spectra_files, database_file, output_dir, output_count)
        os.system(index_cmd)
        os.system(search_cmd)
        output_count += 1
