import numpy as np
import pandas as pd


def print_matrix(matrix, labels):
    """Print the current distance matrix with labels."""
    df = pd.DataFrame(matrix, index=labels, columns=labels)
    print(df)
    print()


def upgma(distance_matrix, labels):
    """
    Perform UPGMA on a half-matrix-form distance matrix.

    :param distance_matrix: 2D numpy array (half matrix representation)
    :param labels: List of labels for the rows/columns
    """
    clusters = {i: [label] for i, label in enumerate(labels)}
    n = len(labels)

    # Convert half-matrix to full matrix for easier computation
    full_matrix = np.zeros((n, n))
    full_matrix[np.triu_indices(n, 1)] = distance_matrix
    full_matrix += full_matrix.T  # Symmetric

    print("Initial Distance Matrix:")
    print_matrix(full_matrix, labels)

    while len(clusters) > 1:
        # Find the smallest non-diagonal value
        i, j = np.unravel_index(np.argmin(full_matrix + np.eye(n) * np.inf), full_matrix.shape)
        if i > j:
            i, j = j, i  # Ensure i < j for consistency

        # Print clusters being merged
        print(f"Merging clusters: {clusters[i]} and {clusters[j]}")

        # Merge clusters
        new_cluster = clusters[i] + clusters[j]
        del clusters[j], clusters[i]
        new_index = len(clusters)
        clusters[new_index] = new_cluster

        # Update the distance matrix
        new_distances = []
        for k in clusters.keys():
            if k != new_index:
                avg_distance = (full_matrix[i, k] + full_matrix[j, k]) / 2
                new_distances.append(avg_distance)

        # Update full matrix
        new_size = len(clusters)
        new_full_matrix = np.zeros((new_size, new_size))
        keys = list(clusters.keys())

        for p, key_p in enumerate(keys[:-1]):
            for q, key_q in enumerate(keys[:-1]):
                new_full_matrix[p, q] = full_matrix[key_p, key_q]

        new_full_matrix[:-1, -1] = new_distances
        new_full_matrix[-1, :-1] = new_distances
        new_full_matrix[-1, -1] = 0  # Self-distance

        # Prepare for next iteration
        full_matrix = new_full_matrix
        n = new_size

        # Print updated distance matrix
        new_labels = [','.join(clusters[k]) for k in keys]
        print("Updated Distance Matrix:")
        print_matrix(full_matrix, new_labels)

    print("Final Cluster:", clusters)


# Example Usage
distance_matrix = np.array([
    7, 10, 14,  # Distances from cluster 1 to clusters 2, 3, 4
    9, 12,  # Distances from cluster 2 to clusters 3, 4
    13  # Distances from cluster 3 to cluster 4
])
labels = ['A', 'B', 'C', 'D', 'E', 'F', 'G']
upgma(distance_matrix, labels)
