from utils import __file_exists

'''from_csv

DESC:
    read in proteins from csv into memory
PARAMS:
    csv_file: str path to csv file
OPTIONAL:
    has_header: bool if set to true, skip header. Assumed in order name, sequence. Default=False
RETURNS:
    list of dictionaries of form {name: str, sequence: str}
'''
def from_csv(csv_file, has_header=False):
    if not __file_exists(csv_file):
        raise Exception('File {} does not exist'.format(csv_file))
    read_header = False
    prots = []
    with open(csv_file, 'r') as i:
        for line in i:
            if has_header and not read_header:
                read_header = True
                continue
            l = line.split(',')
            prots.append({
                'name': l[0],
                'sequence': l[1]
            })
    return prots
