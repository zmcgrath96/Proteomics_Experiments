from utils import __make_valid_csv_file, __make_valid_dir_string, __make_dir

'''write_iter_of_dicts

DESC:
    write a csv from an iterable of dictionaries with entries as column names
PARAMS:
    iter_of_dicts: an iterable of dictionaries to write
    file_name: name of file to write to
OPTIONAL:
    no_header: bool don't use a header line. Default=False
RETURNS: 
    none
'''
def write_iter_of_dicts(iter_of_dicts, file_name, no_header=False):
    if len(iter_of_dicts) == 0:
        return 

    d = '/'.join(file_name.split('/')[:-1])
    d = __make_valid_dir_string(d)
    __make_dir(d)
    file_name = __make_valid_csv_file(file_name)
    d = None 
    keys = [key for key in iter_of_dicts[0]]
    template = ''
    for _ in range(len(keys)):
        template += '{},'
    template = template[:-1] + '\n'
    with open(file_name, 'w') as o:
        if not no_header:
            o.write(template.format(*keys))
        for d in iter_of_dicts:
            vals = [i for _, i in d.items()]
            o.write(template.format(*vals))     

    return(file_name) 


