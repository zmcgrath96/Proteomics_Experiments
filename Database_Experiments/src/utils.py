import os, gzip, shutil

'''__get_related_files

DESC:
    Given some substring, return all files with that substring
PARAMS:
    files: list of strings of file names
    sub: string substring to find in the files
OPTIONAL: 
    not_sub: string a substring that if found in file name don't add
'''
def __get_related_files(files, sub, not_sub=None):
    if not_sub is not None and not_sub != '':
        return [x for x in files if sub in x and not_sub not in x]
    return [x for x in files if sub in x]

'''__make_valid_dir_string

DESC:
    add / to end of director string if it doesn't have it alread
PARAMS:
    dir_path: string path to directory
RETURNS:
    dir path with / at end
'''
def __make_valid_dir_string(dir_path):
    return dir_path + '/' if dir_path[-1] != '/' else dir_path

'''__make_dir

DESC:
    Given a path to directory, check if it exists and if not create it
PARAMS:
    dir_path: string path to a directory to make or check
RETURNS:
    None
'''
def __make_dir(dir_path):
    dir_path = __make_valid_dir_string(dir_path)
    if not os.path.exists(dir_path): 
        os.makedirs(dir_path)

'''__make_valid_text_file

DESC:
    make a string into the name for a text file and make sure directory exists
PARAMS:
    file_name: string a name of the file to save
RETURNS:
    file name with .txt after
'''
def __make_valid_text_file(file_name):
    file_name = file_name + '.txt' if '.txt' not in file_name else file_name
    return file_name

'''__make_valid_json_file

DESC:
    make a string into the name for a text file and make sure the directory exists
PARAMS:
    file_name: string of the file to save
RETURNS:
    file name with .json after it
'''
def __make_valid_json_file(file_name):
    file_name = file_name + '.json' if '.json' not in file_name else file_name
    return file_name

'''__make_valid_csv_file

DESC:
    make a string into the name for a text file 
PARAMS:
    file_name: string of the file to save
RETURNS:
    file name with .csv after it
'''
def __make_valid_csv_file(file_name):
    file_name = file_name + '.csv' if '.csv' not in file_name else file_name
    return file_name

'''__make_valid_fasta_file

DESC:
    make a string into the name for a text file 
PARAMS:
    file_name: string of the file to save
RETURNS:
    file name with .fasta after it
'''
def __make_valid_fasta_file(file_name):
    file_name = file_name + '.fasta' if '.fasta' not in file_name else file_name
    return file_name

'''__file_exists

DESC:
    find out if a file exists
PARAMS:
    file_name: string name of the file to check for
RETURNS:
    bool true if file exists false otherwise
'''
def __file_exists(file_name):
    return os.path.isfile(file_name)

'''__gzip

DESC:
    zip up a file
PARAMS: 
    file_name: str path to the file name to compress
OPTIONAL:
    delete_old: bool delete the uncompressed file. Default=True
RETURNS:
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
PARAMS:
    compressed_file_name: str name of the compressed file to unzip
OPTIONAL:
    delete_old: bool delete the compressed file. Default=True
RETURNS:
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
PARAMS:
    file_name: str path to file in question
RETURNS:
    bool True if file is compressed else False
'''
def __is_gzipped(file_name):
    return '.gz' == file_name[-3:]

'''__gzip_dir

DESC:
    compress a directory with gzip
PARAMS:
    d: str path to directory
OPTIONAL:
    delete_old: bool delete the unziped directory. Default=True
RETURNS:
    path to the new zipped folder
'''
def __gzip_dir(d, delete_old=True):
    root = '/'.join(d.split('/')[:-1])
    shutil.make_archive(d, 'zip', root)
    delete_old and shutil.rmtree(d)
    return d + '.zip'