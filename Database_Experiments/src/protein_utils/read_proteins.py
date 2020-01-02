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

'''__is_dup

DESC:
    see if the new entry is a duplicate
PARAMS: 
    tracker: dictionary {name: seq} to determine keep track
    entry: dictionary {name: str, seq: str} the new entry
RETURNS:
    True if the new entry has been added before, false otherwise
'''
def __is_dup(tracker, entry):
    if entry['name'] in tracker:
        print('repeat name found: {}'.format(entry['name']))
    return entry['name'] in tracker and entry['seq'] == tracker[entry['name']]

'''from_fasta

DESC:
    read proteins into memory from fasta file
PARAMS:
    fasta_file: str path to fasta file
RETURNS:
    list of dictionaries of form {'name': str, 'sequence': str, 'identifier': str}
    list of dictionaries of duplicates of the same form as above
'''
def from_fasta(fasta_file, save_dir='./'):
    if not __file_exists(fasta_file):
        raise Exception('File {} does not exist'.format(fasta_file))
    prots = []
    dups = []
    tracker = {}
    print('Loading proteins...')
    with open(fasta_file, 'r') as i:
        name = None 
        seq = '' 
        identifier = ''
        for line in i:
            if '>' in line: #name line

                # add the last thing to the list
                if not ((name is None or name == '') and (seq is None or seq == '')):
                    if __is_dup(tracker, {'name': name, 'seq': seq}):
                        dups.append({
                            'name': name,
                            'sequence': seq,
                            'identifier': identifier
                        })

                    tracker[name] = seq
                    prots.append({
                        'name': name,
                        'sequence': seq,
                        'identifier': identifier
                    })

                seq = '' 
                name = str(str(line.split('|')[2]).split(' ')[0]).replace('\n', '')
                identifier = str(line.split('|')[1])
            else:
                seq += line.replace('\n', '')
        # add the last one
        if __is_dup(tracker, {'name': name, 'seq': seq}):
            dups.append({
                'name': name,
                'sequence': seq,
                'identifier': identifier
            })
        else:
            prots.append({
                'name': name,
                'sequence': seq,
                'identifier': identifier
            })
    print('Done.')
    return prots, dups
