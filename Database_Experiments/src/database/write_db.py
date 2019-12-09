from utils import __make_valid_dir_string, __make_dir

def write_fasta(output_path, output_name, sequences, prefix):
    if '.fasta' not in output_name:
        output_name += '.fasta'

    output_path = __make_valid_dir_string(output_path)
    __make_dir(output_path)
    output_file = output_path + output_name
    count = 0
    with open(output_file, 'w') as o:
        for sequence in sequences:
            o.write('>' + prefix + str(count) + '\n')
            count += 1
            o.write(sequence + '\n')

    return output_file