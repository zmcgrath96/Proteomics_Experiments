import pyopenms
import os
from subprocess import call
from utils import __make_dir, __make_valid_dir_string, __is_gzipped, __gunzip, __gzip

crux_to_rm = ['tide-search.decoy.txt', 'tide-search.log.txt', 'tide-search.params.txt']

def __parse_spectrum_name(spec_name):
    return str(spec_name.split('/')[-1]).lower().replace('.mzml', '')

def __parse_db_name(db_name):
    return  str(db_name.split('/')[-1]).replace('.fasta', '')

'''score_peptides

DESC:
    use the crux tool to score spectra against databases
PARAMS:
    spectra_files: list of str paths to all the spectra (.mzML) files
    database_files: list of str paths to all the database (.fasta) files
    path_to_crux_cmd: str path to the executable for crux
    output_dir: str path to the directory to save files
OPTIONAL:
    compress: bool compress the output result. Default=True
RETURNS:
    list of str of output files
'''
def score_peptides(spectra_files, database_files, path_to_crux_cmd, output_dir, compress=True):
    output_dir = __make_valid_dir_string(output_dir) + 'search_output/'
    __make_dir(output_dir)

    is_compressed = __is_gzipped(spectra_files[0])

    output_count = 0
    output_files = []
    num_dbs = len(database_files)
    num_specs = len(spectra_files)

    for i, spec_file in enumerate(spectra_files):
        spec_file = spec_file if not is_compressed else __gunzip(spec_file)
        for j, database_file in enumerate(database_files):
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
            output_files.append(o)

            # saving space, so should remove the extra stuff
            if is_compressed:
                for rm in crux_to_rm:
                    os.remove(this_output_dir + '/' + rm)

        # if the files were compressed, we were trying to save disk space so just remove the mzml files
        is_compressed and os.remove(spec_file)

    return output_files