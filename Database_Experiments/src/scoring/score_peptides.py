import pyopenms
import os
from subprocess import call

def __parse_spectrum_name(spec_name):
    end = spec_name.split('/')[-1]
    spl = end.split('_')[1:]
    return '_'.join(spl[:-1] + [spl[-1].replace('.mzml', '')])  

def __parse_db_name(db_name):
    return  str(db_name.split('/')[-1]).replace('.fasta', '')

def score_peptides(spectra_files, database_files, path_to_crux_cmd, output_dir):
    output_dir = output_dir + '/' if output_dir[-1] != '/' else output_dir
    # str_spectra_files = ' '.join(spectra_files)
    output_count = 0
    output_files = []
    num_dbs = len(database_files)
    num_specs = len(spectra_files)
    #index_prefix = 'index_'
    for i, database_file in enumerate(database_files):
        this_db_name = __parse_db_name(database_file)
        for j, spec_file in enumerate(spectra_files):
            print('On database file {}/{}[{}%]\tOn spectrum {}/{}[{}%]\r'.format(i+1, num_dbs, int(((i+1)/num_dbs) * 100), j+1, num_specs, int(((j+1)/num_specs)*100)), end="")
            this_output_dir = output_dir + '{}_vs_{}'.format(__parse_spectrum_name(spec_file), this_db_name)
            #curr_index = index_prefix + str(output_count)
            #index_cmd = [path_to_crux_cmd, 'tide-index', database_file, curr_index, '--min-length', '3', '--overwrite', 'T', '--min-mass', '100'] #'{} tide-index {} {} --min-length 3 --overwrite T'.format(path_to_crux_cmd, database_file, curr_index)
            search_cmd = [
                path_to_crux_cmd, 
                'tide-search', 
                spec_file, 
                database_file, 
                '--min-length', '3', 
                '--min-mass', '100', 
                '--output-dir', this_output_dir, 
                '--overwrite', 'T', 
                '--min-peaks', '2', 
                '--precursor-window', '1000000000',
                '--enzyme', 'no-enzyme', 
                '--verbosity', '0'
                ] #'{} tide-search {} {} --output-dir {}_{} --overwrite T --min-peaks 2'.format(path_to_crux_cmd, str_spectra_files, database_file, output_dir, output_count)
            #call(index_cmd)
            call(search_cmd)
            output_count += 1
            output_files.append(this_output_dir + '/tide-search.target.txt')

    return output_files