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

'''from_fasta

DESC:
    read proteins into memory from fasta file
PARAMS:
    fasta_file: str path to fasta file
RETURNS:
    list of dictionaries of form {'name': str, 'sequence': str}
'''
def from_fasta(fasta_file):
    if not __file_exists(fasta_file):
        raise Exception('File {} does not exist'.format(fasta_file))
    prots = []
    with open(fasta_file, 'r') as i:
        name = None 
        seq = '' 
        for line in i:
            if '>' in line: #name line
                name is not None and prots.append({
                    'name': name,
                    'sequence': seq
                })
                seq = '' 
                name = str(line.replace('>', '').split('|')[0]).replace('\n', '')
            else:
                seq += line.replace('\n', '')
        # add the last one
        prots.append({
            'name': name,
            'sequence': seq
        })
    return prots
