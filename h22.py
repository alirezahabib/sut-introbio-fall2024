import numpy as np


def compute_q_matrix(matrix, total_distances):
    n = matrix.shape[0]
    q_matrix = np.zeros((n, n))

    for i in range(n):
        for j in range(i + 1, n):
            q_value = (n - 2) * matrix[i, j] - total_distances[i] - total_distances[j]
            q_matrix[i, j] = q_matrix[j, i] = q_value

    np.fill_diagonal(q_matrix, float('inf'))
    return q_matrix


def compute_branch_lengths(matrix, i, j, total_distances):
    n = matrix.shape[0]
    dist_ij = matrix[i, j]

    branch_i = 0.5 * (dist_ij + (total_distances[i] - total_distances[j]) / (n - 2))
    branch_j = dist_ij - branch_i

    return max(0, branch_i), max(0, branch_j)


class NeighborJoiningTree:
    def __init__(self, distance_matrix):
        self.n = len(distance_matrix)
        self.distances = np.array(distance_matrix, dtype=float)
        self.tree = {}
        self.current_vertex = self.n

    def neighbor_joining(self):
        matrix = self.distances.copy()
        active_nodes = list(range(self.n))

        while len(active_nodes) > 2:
            # Row totals and Q-matrix
            total_distances = np.sum(matrix, axis=1)
            q_matrix = compute_q_matrix(matrix, total_distances)

            # Minimum Q-matrix indices
            i, j = np.unravel_index(np.argmin(q_matrix), q_matrix.shape)
            i, j = active_nodes[i], active_nodes[j]

            branch_length_i, branch_length_j = compute_branch_lengths(
                matrix,
                active_nodes.index(i),
                active_nodes.index(j),
                total_distances
            )

            # Connecting i and j
            self.tree[i] = self.tree.get(i, []) + [(self.current_vertex, round(branch_length_i, 3))]
            self.tree[j] = self.tree.get(j, []) + [(self.current_vertex, round(branch_length_j, 3))]
            self.tree[self.current_vertex] = [
                (i, round(branch_length_i, 3)),
                (j, round(branch_length_j, 3))
            ]

            # New distance matrix row for the merged vertex
            new_distances = 0.5 * (matrix[active_nodes.index(i)] + matrix[active_nodes.index(j)])

            # Remove the merged nodes
            active_indices = list(range(len(active_nodes)))
            active_indices.remove(active_nodes.index(i))
            active_indices.remove(active_nodes.index(j))

            matrix = matrix[active_indices][:, active_indices]

            # Append new distances
            matrix = np.column_stack([matrix, new_distances[active_indices]])
            matrix = np.row_stack([matrix, np.append(new_distances[active_indices], 0)])

            # Update active nodes
            active_nodes = [node for node in active_nodes if node not in [i, j]]
            active_nodes.append(self.current_vertex)

            self.current_vertex += 1

        # Final two nodes
        if len(active_nodes) == 2:
            i, j = active_nodes
            final_dist = matrix[0, 1]
            self.tree[i] = self.tree.get(i, []) + [(j, round(final_dist, 3))]
            self.tree[j] = self.tree.get(j, []) + [(i, round(final_dist, 3))]

        return self.tree


def print_adjacency_list(tree):
    for node, neighbors in sorted(tree.items()):
        for neighbor, weight in sorted(neighbors):
            print(f'{node}->{neighbor}:{weight:.3f}')


n = int(input())
matrix = []
for _ in range(n):
    matrix.append(list(map(float, input().split())))

nj_tree = NeighborJoiningTree(matrix)
tree = nj_tree.neighbor_joining()
print_adjacency_list(tree)
