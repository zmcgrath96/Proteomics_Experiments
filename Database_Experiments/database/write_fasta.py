def write(output_path, output_name, sequences, prefix):
    if '.fasta' not in output_name:
        output_name += '.fasta'

    output_file = output_path + '/' + output_name
    count = 0
    with open(output_file, 'w') as o:
        for sequence in sequences:
            o.write('>' + prefix + str(count) + '\n')
            count += 1
            o.write(sequence + '\n')