def k_mers(original_sequence, window_size):
    if window_size > len(original_sequence):
        return [original_sequence]
        
    w_1 = window_size-1
    return [original_sequence[i:i+window_size] for i in range(len(original_sequence)-w_1)]