import os, gzip, shutil, copy

'''__get_related_files

DESC:
    Given some substring, return all files with that substring
Inputs:
    files: list of strings of file names
    sub: string substring to find in the files
kwargs: 
    not_sub: string a substring that if found in file name don't add
'''
def __get_related_files(files, sub, not_sub=None):
    if not_sub is not None and not_sub != '':
        return [x for x in files if sub in x and not_sub not in x]
    return [x for x in files if sub in x]

'''__make_valid_dir_string

DESC:
    add / to end of director string if it doesn't have it alread
Inputs:
    dir_path: string path to directory
Outputs:
    dir path with / at end
'''
def __make_valid_dir_string(dir_path):
    return dir_path + '/' if dir_path[-1] != '/' else dir_path

'''__make_dir

DESC:
    Given a path to directory, check if it exists and if not create it
Inputs:
    dir_path: string path to a directory to make or check
Outputs:
    None
'''
def __make_dir(dir_path):
    dir_path = __make_valid_dir_string(dir_path)
    if not os.path.exists(dir_path): 
        os.makedirs(dir_path)

'''__make_valid_text_file

DESC:
    make a string into the name for a text file and make sure directory exists
Inputs:
    file_name: string a name of the file to save
Outputs:
    file name with .txt after
'''
def __make_valid_text_file(file_name):
    file_name = file_name + '.txt' if '.txt' not in file_name else file_name
    return file_name

'''__make_valid_json_file

DESC:
    make a string into the name for a text file and make sure the directory exists
Inputs:
    file_name: string of the file to save
Outputs:
    file name with .json after it
'''
def __make_valid_json_file(file_name):
    file_name = file_name + '.json' if '.json' not in file_name else file_name
    return file_name

'''__make_valid_csv_file

DESC:
    make a string into the name for a text file 
Inputs:
    file_name: string of the file to save
Outputs:
    file name with .csv after it
'''
def __make_valid_csv_file(file_name):
    file_name = file_name + '.csv' if '.csv' not in file_name else file_name
    return file_name

'''__make_valid_fasta_file

DESC:
    make a string into the name for a text file 
Inputs:
    file_name: string of the file to save
Outputs:
    file name with .fasta after it
'''
def __make_valid_fasta_file(file_name):
    file_name = file_name + '.fasta' if '.fasta' not in file_name else file_name
    return file_name

'''__file_exists

DESC:
    find out if a file exists
Inputs:
    file_name: string name of the file to check for
Outputs:
    bool true if file exists false otherwise
'''
def __file_exists(file_name):
    return os.path.isfile(file_name)

'''__gzip

DESC:
    zip up a file
Inputs: 
    file_name: str path to the file name to compress
kwargs:
    delete_old: bool delete the uncompressed file. Default=True
Outputs:
    str name of the new compressed file
'''
def __gzip(file_name, delete_old=True):
    compressed_file_name = file_name + '.gz'
    with open(file_name, 'rb') as f_in:
        with gzip.open(compressed_file_name, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    delete_old and os.remove(file_name)
    return compressed_file_name

'''__gunzip

DESC:   
    unzip a file
Inputs:
    compressed_file_name: str name of the compressed file to unzip
kwargs:
    delete_old: bool delete the compressed file. Default=True
Outputs:
    str name of the file unziped
'''  
def __gunzip(compressed_file_name, delete_old=True):
    file_name = compressed_file_name if '.gz' not in compressed_file_name else compressed_file_name.replace('.gz', '')
    with gzip.open(compressed_file_name, 'rb') as f_in:
        with open(file_name, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    delete_old and os.remove(compressed_file_name)
    return file_name

'''__is_gzipped

DESC:
    determines if a file has been gzipped
Inputs:
    file_name: str path to file in question
Outputs:
    bool True if file is compressed else False
'''
def __is_gzipped(file_name):
    return '.gz' == file_name[-3:]

'''__gzip_dir

DESC:
    compress a directory with gzip
Inputs:
    d: str path to directory
kwargs:
    delete_old: bool delete the unziped directory. Default=True
Outputs:
    path to the new zipped folder
'''
def __gzip_dir(d, delete_old=True):
    root = '/'.join(d.split('/')[:-1])
    shutil.make_archive(d, 'zip', root)
    delete_old and shutil.rmtree(d)
    return d + '.zip'

'''__is_json

DESC:
    determine if a file is a json file based purely on name
Inputs:
    file: file to determine if its a json file
Outputs:
    bool True if it is a json file False otherwise
'''
def __is_json(file):
    return True if '.json' in file else False

'''__is_fasta

DESC:
    determine if a file is a fasta file based purely on name
Inputs:
    file: file to determine if its a fasta file
Outputs:
    bool True if it is a fasta file False otherwise
'''
def __is_fasta(file):
    return True if '.fasta' in file else False

def __split_exp_by_ion(exp: dict, ion: str) -> dict:
    '''
    Make a copy of the experiment split by ion information

    Inputs:
        exp:    dictionary with expeiment summary information
        ion:    string ion type. Possible types are {'b', 'y'}
    Outpus:
        dict same structure as experiment just without other ion information
    '''
    ionized_exp = {}
    ionized_exp['experiment_info'] = copy.deepcopy(exp['experiment_info'])
    ionized_exp['experiment'] = {}
    for pep_name, pep in exp['experiment'].items():
        ionized_exp['experiment'][pep_name] = {}
        for prot_name, prot in pep.items():
            if prot_name == 'analysis':
                ionized_exp['experiment'][pep_name][prot_name] = copy.deepcopy(prot)
                ionized_exp['experiment'][pep_name][prot_name]['ranks']['ranks'] = ionized_exp['experiment'][pep_name][prot_name]['ranks']['ranks'][ion]
                continue
            ionized_exp['experiment'][pep_name][prot_name] = {}
            for k, scores in prot.items():
                ionized_exp['experiment'][pep_name][prot_name][k] = scores[ion] 

    return ionized_exp