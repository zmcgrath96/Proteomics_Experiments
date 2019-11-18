def gen_sequences(original_sequence, window_size):
    seqs = []
    if window_size > len(original_sequence):
        seqs = [original_sequence]
        return seqs
    for i in range(len(original_sequence) - window_size + 1):
        seqs.append(original_sequence[i:i+window_size])

    return seqs