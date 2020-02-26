from utils import __make_valid_fasta_file, __file_exists

'''write

DESC:
    write a fasta file
Inputs:
    output_name: str name of file to write to 
    sequences: list of dictionaries of form {'name': str, 'sequence': str}
Outputs:
    name of the output file written to
'''
def write(output_name, sequences):
    output_name = __make_valid_fasta_file(output_name)
    with open(output_name, 'w') as o:
        for i, seq in enumerate(sequences):
            o.write('>sp|{}|{}\n{}\n'.format('id{}'.format(i), seq['name'], seq['sequence']))
    return output_name

'''read

DESC:
    read proteins into memory from fasta file
Inputs:
    fasta_file: str path to fasta file
Outputs:
    list of dictionaries of form {'name': str, 'sequence': str, 'identifier': str}
'''
def read(fasta_file):
    if not __file_exists(fasta_file):
        raise Exception('File {} does not exist'.format(fasta_file))
    prots = []
    with open(fasta_file, 'r') as i:
        name = None 
        seq = '' 
        identifier = ''
        for line in i:
            if '>' in line: #name line

                # add the last thing to the list
                if not ((name is None or name == '') and (seq is None or seq == '')):
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
        prots.append({
            'name': name,
            'sequence': seq,
            'identifier': identifier
        })
    return prots