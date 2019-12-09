import os 

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
    d = '/'.join(file_name.split('/')[:-1])
    __make_dir(d)
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
    d = '/'.join(file_name.split('/')[:-1])
    __make_dir(d)
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