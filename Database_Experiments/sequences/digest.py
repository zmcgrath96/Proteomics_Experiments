from random import randint, choice, random

'''tryptic

DESC:
    Generate peptides with tryptic digest from sequences with missed cleavages at a probability
'''
def tryptic(sequences, number_digests, miss_prob=0, save_dir='./', save_name='digestion.tsv'):
    to_digest = []
    if number_digests <= len(sequences):
        for _ in range(number_digests):
            to_digest.append(choice(sequences))
    else:
        to_digest += sequences
        for _ in range(number_digests - len(sequences)):
            to_digest.append(choice(sequences))

    peptides = []
    o = open(save_dir + save_name, 'w')
    form = '{}\t{}\t{}\t{}\n'
    o.write(form.format('peptide', 'parent-name', 'parent-sequence', 'start-location'))

    for digest in to_digest:
        seq = digest['sequence']
        name = digest['name']
        this_pep = ''
        # pick a start point thats not the start or the end
        start = randint(0, len(seq)-1)
        while start == 0 or start == len(seq) - 1:
            start = randint(0, len(seq) - 1)
        end = start

        # if start is at a seq point, go left or right
        if seq[start] == 'R' or seq[start] == 'K':
            # go left
            if random() < 0.5:
                while start > 0:
                    start -= 1
                    if seq[start] == 'K' or seq[start] == 'R':
                        if random() < miss_prob:
                            continue
                        break 
                this_pep = seq[start+1:end+1]
            # go right
            else:
                while end < len(seq) - 1:
                    end += 1
                    if seq[end] == 'K' or seq[end] == 'R':
                        if random() < miss_prob:
                            continue
                        break 
                this_pep = seq[start+1:end+1]

        # go both ways
        else:
            while start > 0:
                if seq[start] == 'K' or seq[start] == 'R':
                    if random() < miss_prob:
                        continue
                    break 
                start -= 1
            while end < len(seq) - 1:
                if seq[end] == 'K' or seq[end] == 'R':
                    if random() < miss_prob:
                        continue
                    break 
                end += 1
            this_pep = seq[start+1:end+1]

        peptides.append(this_pep)
        o.write(form.format(this_pep, name, seq, start))

    return peptides
