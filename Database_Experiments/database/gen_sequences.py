def gen_sequences(original_sequence, window_size):
    seqs = []
    for i in range(len(original_sequence) - window_size + 1):
        seqs.append(original_sequence[i:i+window_size])

    return seqs