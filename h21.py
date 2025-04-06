"""
Write Python code for this problem. Researchers are analyzing a number of genes to solve a problem. Unfortunately, the exact sequences of these genes have been lost. However, all the k-mers (substrings of length k) of these sequences (including duplicate k-mers) have been preserved. The researchers need to determine how many times a query string appears as a substring in one of the possible input genes. To ensure accuracy, the researchers want to calculate the maximum possible occurrences of the query string among all possible sequence sets that match the given k-mer data.

# **Input** 

The first line of input contains two integers, n and k, where n is the number of k-mers and k is the length of each k-mer. The next n lines each contain one k-mer. After that, an integer q indicates the number of queries, followed by q lines, each containing one query string.

# **Output** 

For each query, output the maximum possible occurrences of the query string as a substring among all valid gene sequences that match the k-mer input. If it is not possible to determine an upper bound on the maximum occurrences, output -1.



In essence, for each query, you must determine the maximum number of occurrences among all valid sequences consistent with the given k-mer set.
"""

# #%%
# #your code
import numpy as np
from collections import Counter

def max_occurrences(kmers, query):
    n = len(kmers)
    k = len(kmers[0])
    query_count = Counter(query)
    kmers_count = [Counter(kmer) for kmer in kmers]
    max_occurrences = 0

    for i in range(n):
        valid = True
        for char, count in query_count.items():
            if kmers_count[i][char] < count:
                valid = False
                break

        if valid:
            occurrences = 1
            for char, count in query_count.items():
                occurrences *= kmers_count[i][char] // count

            max_occurrences = max(max_occurrences, occurrences)

    return max_occurrences if max_occurrences else -1

# Sample Input
n, k = map(int, input().split())

kmers = [input() for _ in range(n)]
q = int(input())
queries = [input() for _ in range(q)]

# Output
for query in queries:
    print(max_occurrences(kmers, query))


