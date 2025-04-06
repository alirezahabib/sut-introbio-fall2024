from numpy import zeros

inf = 1e6

BLOSUM62 = [
    [4, 0, -2, -1, -2, 0, -2, -1, -1, -1, -1, -2, -1, -1, -1, 1, 0, 0, -3, -2],
    [0, 9, -3, -4, -2, -3, -3, -1, -3, -1, -1, -3, -3, -3, -3, -1, -1, -1, -2, -2],
    [-2, -3, 6, 2, -3, -1, -1, -3, -1, -4, -3, 1, -1, 0, -2, 0, -1, -3, -4, -3],
    [-1, -4, 2, 5, -3, -2, 0, -3, 1, -3, -2, 0, -1, 2, 0, 0, -1, -2, -3, -2],
    [-2, -2, -3, -3, 6, -3, -1, 0, -3, 0, 0, -3, -4, -3, -3, -2, -2, -1, 1, 3],
    [0, -3, -1, -2, -3, 6, -2, -4, -2, -4, -3, 0, -2, -2, -2, 0, -2, -3, -2, -3],
    [-2, -3, -1, 0, -1, -2, 8, -3, -1, -3, -2, 1, -2, 0, 0, -1, -2, -3, -2, 2],
    [-1, -1, -3, -3, 0, -4, -3, 4, -3, 2, 1, -3, -3, -3, -3, -2, -1, 3, -3, -1],
    [-1, -3, -1, 1, -3, -2, -1, -3, 5, -2, -1, 0, -1, 1, 2, 0, -1, -2, -3, -2],
    [-1, -1, -4, -3, 0, -4, -3, 2, -2, 4, 2, -3, -3, -2, -2, -2, -1, 1, -2, -1],
    [-1, -1, -3, -2, 0, -3, -2, 1, -1, 2, 5, -2, -2, 0, -1, -1, -1, 1, -1, -1],
    [-2, -3, 1, 0, -3, 0, 1, -3, 0, -3, -2, 6, -2, 0, 0, 1, 0, -3, -4, -2],
    [-1, -3, -1, -1, -4, -2, -2, -3, -1, -3, -2, -2, 7, -1, -2, -1, -1, -2, -4, -3],
    [-1, -3, 0, 2, -3, -2, 0, -3, 1, -2, 0, 0, -1, 5, 1, 0, -1, -2, -2, -1],
    [-1, -3, -2, 0, -3, -2, 0, -3, 2, -2, -1, 0, -2, 1, 5, -1, -1, -3, -3, -2],
    [1, -1, 0, 0, -2, 0, -1, -2, 0, -2, -1, 1, -1, 0, -1, 4, 1, -2, -3, -2],
    [0, -1, -1, -1, -2, -2, -2, -1, -1, -1, -1, 0, -1, -1, -1, 1, 5, 0, -2, -2],
    [0, -1, -3, -2, -1, -3, -3, 3, -2, 1, 1, -3, -2, -2, -3, -2, 0, 4, -3, -1],
    [-3, -2, -4, -3, 1, -2, -2, -3, -3, -2, -1, -4, -4, -2, -3, -3, -2, -3, 11, 2],
    [-2, -2, -3, -2, 3, -3, 2, -1, -2, -1, -1, -2, -3, -1, -2, -2, -2, -1, 2, 7]
]

mapping = {a: i for i, a in enumerate('ACDEFGHIKLMNPQRSTVWY')}


def global_align(s1, s2, scores, gap_pen, ext_pen):
    # Score
    M = zeros((len(s1) + 1, len(s2) + 1), dtype=int)  # match/mismatch
    X = zeros((len(s1) + 1, len(s2) + 1), dtype=int)  # X gap
    Y = zeros((len(s1) + 1, len(s2) + 1), dtype=int)  # Y gap

    # Traceback
    M_trace = zeros((len(s1) + 1, len(s2) + 1), dtype=int)
    X_trace = zeros((len(s1) + 1, len(s2) + 1), dtype=int)
    Y_trace = zeros((len(s1) + 1, len(s2) + 1), dtype=int)

    # First row/column
    M[1:len(s1) + 1, 0] = [gap_pen + ext_pen * (i - 1) for i in range(1, len(s1) + 1)]
    M[0, 1:len(s2) + 1] = [gap_pen + ext_pen * (i - 1) for i in range(1, len(s2) + 1)]

    X[1:len(s1) + 1, 0] = Y[1:len(s1) + 1, 0] = X[0, 1:len(s2) + 1] = Y[0, 1:len(s2) + 1] = -inf

    # Filling
    for i in range(1, len(s1) + 1):
        for j in range(1, len(s2) + 1):
            X[i][j], X_trace[i][j] = max((M[i - 1][j] + gap_pen, 0), (X[i - 1][j] + ext_pen, 1))
            Y[i][j], Y_trace[i][j] = max((M[i][j - 1] + gap_pen, 0), (Y[i][j - 1] + ext_pen, 1))
            M[i][j], M_trace[i][j] = max((M[i - 1][j - 1] + scores[mapping[s1[i - 1]]][mapping[s2[j - 1]]], 0),
                                         (X[i][j], 1), (Y[i][j], 2))

    s1p, s2p = s1, s2
    best_score, trace = max((X[-1][-1], 0), (Y[-1][-1], 1), (M[-1][-1], 2))
    i, j = len(s1), len(s2)

    # Traceback to align
    while i > 0 and j > 0:
        if trace == 0:
            if X_trace[i][j] == 0:
                trace = 2
            i -= 1
            s2p = s2p[:j] + '-' + s2p[j:]
        elif trace == 1:
            if Y_trace[i][j] == 0:
                trace = 2
            j -= 1
            s1p = s1p[:i] + '-' + s1p[i:]
        else:
            if M_trace[i][j] == 1:
                trace = 0
            elif M_trace[i][j] == 2:
                trace = 1
            else:
                i -= 1
                j -= 1

    s2p = '-' * i + s2p
    s1p = '-' * j + s1p

    return str(best_score), s1p, s2p


def parse_fasta(fasta_file):
    with open(fasta_file, 'r') as file:
        return [''.join(entry.splitlines()[1:]) for entry in file.read().split('>') if entry]


dna1, dna2 = parse_fasta('rosalind_gaff.txt')
print('\n'.join(global_align(dna1, dna2, BLOSUM62, -11, -1)))
