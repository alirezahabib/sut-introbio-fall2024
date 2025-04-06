from numpy import zeros

PAM250 = [
	[ 2, -2,  0,  0, -3,  1, -1, -1, -1, -2, -1,  0,  1,  0, -2,  1,  1,  0, -6, -3],
	[-2, 12, -5, -5, -4, -3, -3, -2, -5, -6, -5, -4, -3, -5, -4,  0, -2, -2, -8,  0],
	[ 0, -5,  4,  3, -6,  1,  1, -2,  0, -4, -3,  2, -1,  2, -1,  0,  0, -2, -7, -4],
	[ 0, -5,  3,  4, -5,  0,  1, -2,  0, -3, -2,  1, -1,  2, -1,  0,  0, -2, -7, -4],
	[-3, -4, -6, -5,  9, -5, -2,  1, -5,  2,  0, -3, -5, -5, -4, -3, -3, -1,  0,  7],
	[ 1, -3,  1,  0, -5,  5, -2, -3, -2, -4, -3,  0,  0, -1, -3,  1,  0, -1, -7, -5],
	[-1, -3,  1,  1, -2, -2,  6, -2,  0, -2, -2,  2,  0,  3,  2, -1, -1, -2, -3,  0],
	[-1, -2, -2, -2,  1, -3, -2,  5, -2,  2,  2, -2, -2, -2, -2, -1,  0,  4, -5, -1],
	[-1, -5,  0,  0, -5, -2,  0, -2,  5, -3,  0,  1, -1,  1,  3,  0,  0, -2, -3, -4],
	[-2, -6, -4, -3,  2, -4, -2,  2, -3,  6,  4, -3, -3, -2, -3, -3, -2,  2, -2, -1],
	[-1, -5, -3, -2,  0, -3, -2,  2,  0,  4,  6, -2, -2, -1,  0, -2, -1,  2, -4, -2],
	[ 0, -4,  2,  1, -3,  0,  2, -2,  1, -3, -2,  2,  0,  1,  0,  1,  0, -2, -4, -2],
	[ 1, -3, -1, -1, -5,  0,  0, -2, -1, -3, -2,  0,  6,  0,  0,  1,  0, -1, -6, -5],
	[ 0, -5,  2,  2, -5, -1,  3, -2,  1, -2, -1,  1,  0,  4,  1, -1, -1, -2, -5, -4],
	[-2, -4, -1, -1, -4, -3,  2, -2,  3, -3,  0,  0,  0,  1,  6,  0, -1, -2,  2, -4],
	[ 1,  0,  0,  0, -3,  1, -1, -1,  0, -3, -2,  1,  1, -1,  0,  2,  1, -1, -2, -3],
	[ 1, -2,  0,  0, -3,  0, -1,  0,  0, -2, -1,  0,  0, -1, -1,  1,  3,  0, -5, -3],
	[ 0, -2, -2, -2, -1, -1, -2,  4, -2,  2,  2, -2, -1, -2, -2, -1,  0,  4, -6, -2],
	[-6, -8, -7, -7,  0, -7, -3, -5, -3, -2, -4, -4, -6, -5,  2, -2, -5, -6, 17,  0],
	[-3,  0, -4, -4,  7, -5,  0, -1, -4, -1, -2, -2, -5, -4, -4, -3, -3, -2,  0, 10]
]

mapping = {a: i for i, a in enumerate('ACDEFGHIKLMNPQRSTVWY')}


def local_alignment(seq1, seq2, scoring_matrix, gap_pen):
    s = zeros((len(seq1) + 1, len(seq2) + 1), dtype=int)
    b = zeros((len(seq1) + 1, len(seq2) + 1), dtype=int)  # backtrack matrix, 0 = up, 1 = left, 2 = diagonal, 3 = stop

    for i in range(1, len(seq1) + 1):
        for j in range(1, len(seq2) + 1):
            up = s[i - 1][j] - gap_pen  # gap in seq2
            left = s[i][j - 1] - gap_pen  # gap in seq1
            diag = s[i - 1][j - 1] + scoring_matrix[mapping[seq1[i - 1]]][mapping[seq2[j - 1]]]  # match/mismatch
            zero = 0  # no alignment

            s[i][j], b[i][j] = max((up, 0), (left, 1), (diag, 2), (zero, 3))

    max_index = s.argmax()
    i, j = divmod(max_index, s.shape[1])
    best_score = str(s[i][j])

    seq1_aligned, seq2_aligned = seq1[:i], seq2[:j]

    while i != 0 and j != 0 and b[i][j] != 3:
        i, j = (i - 1, j) if b[i][j] == 0 else (i, j - 1) if b[i][j] == 1 else (i - 1, j - 1)

    return best_score, seq1_aligned[i:], seq2_aligned[j:]

def parse_fasta(fasta_file):
    with open(fasta_file, 'r') as file:
        return [''.join(entry.splitlines()[1:]) for entry in file.read().split('>') if entry]


dna1, dna2 = parse_fasta('rosalind_loca.txt')
print('\n'.join(local_alignment(dna1, dna2, PAM250, 5)))
