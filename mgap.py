import numpy as np

inf = 1e6


def max_gap(s1, s2):
    s = np.zeros((len(s1) + 1, len(s2) + 1), dtype=int)
    gap_count = np.zeros((len(s1) + 1, len(s2) + 1), dtype=int)
    s[1:, 0] = -np.arange(1, len(s1) + 1)
    gap_count[1:, 0] = np.arange(1, len(s1) + 1)
    s[0, 1:] = -np.arange(1, len(s2) + 1)
    gap_count[0, 1:] = np.arange(1, len(s2) + 1)

    for i in range(1, len(s1) + 1):
        for j in range(1, len(s2) + 1):
            s[i][j], gap_count[i][j] = max(
                (s[i - 1][j - 1] + (1 if s1[i - 1] == s2[j - 1] else -inf), gap_count[i - 1][j - 1]),
                (s[i - 1][j] - 1, gap_count[i - 1][j] + 1),
                (s[i][j - 1] - 1, gap_count[i][j - 1] + 1))
    return gap_count[-1][-1]


def parse_fasta(fasta_file):
    with open(fasta_file, 'r') as file:
        return [''.join(entry.splitlines()[1:]) for entry in file.read().split('>') if entry]


dna1, dna2 = parse_fasta('rosalind_mgap.txt')
print(max_gap(dna1, dna2))
