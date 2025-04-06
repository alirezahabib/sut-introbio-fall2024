import numpy as np

def overlap_align(s1, s2):
    s = np.zeros((len(s1) + 1, len(s2) + 1), dtype=int)
    trace = np.zeros((len(s1) + 1, len(s2) + 1), dtype=int)

    # Filling
    for i in range(1, len(s1) + 1):
        for j in range(1, len(s2) + 1):
            s[i][j], trace[i][j] = max((s[i - 1][j - 1] + (1 if s1[i - 1] == s2[j - 1] else -2), 0),
                                       (s[i - 1][j] - 2, 1),
                                       (s[i][j - 1] - 2, 2))

    best_row = np.argmax(s[-1])  # Max value index in the last row
    best_col = np.argmax(s[:, -1])  # Max value index in the last column

    if best_row > best_col:
        i, j = len(s1), best_row
    else:
        i, j = best_col, len(s2)

    best_score = s[i][j]
    s1p, s2p = s1[:i], s2[:j]

    # Traceback
    while i > 0 and j > 0:
        if trace[i][j] == 0:
            i -= 1
            j -= 1
        elif trace[i][j] == 1:
            i -= 1
            s2p = s2p[:j] + '-' + s2p[j:]
        else:
            j -= 1
            s1p = s1p[:i] + '-' + s1p[i:]


    s1p = s1p[i:]
    return str(best_score), s1p, s2p


print('\n'.join(overlap_align(input(), input())))
