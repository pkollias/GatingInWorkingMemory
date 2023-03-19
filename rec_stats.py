from rec import TimebinInterval
from versioning import *
from rec_db import *
from sklearn.preprocessing import StandardScaler
from math import degrees

def estimate_cluster_mass_2d(cluster_vals_array, alternative='greater', observed=False):

    nrows, ncols = cluster_vals_array.shape
    if alternative == 'greater':
        island_mask = cluster_vals_array > 0
    elif alternative == 'less':
        island_mask = cluster_vals_array < 0

    def absorb_cluster(row, col, array, mask, inds_list, observed=False):
        nrows, ncols = array.shape
        if row < 2 or col < 2 or row >= 41 or col >= 41:
            return 0
        elif not mask[row][col]:
            return 0
        else:
            if observed:
                inds_list.append((row, col))
            mask[row][col] = False
            neighbor_list = [(row + 1, col), (row , col + 1), (row - 1, col), (row, col - 1)]
            return array[row][col] + sum([absorb_cluster(row_i, col_j, array, mask, inds_list, observed) for row_i, col_j in neighbor_list])

    edge_cluster = 0
    edge_cluster_list = []
    for row, col in product(range(nrows), range(ncols)):
        inds_list = []
        cluster_mass = absorb_cluster(row, col, cluster_vals_array, island_mask, inds_list, observed)
        if (cluster_mass > edge_cluster and alternative == 'greater') or \
           (cluster_mass < edge_cluster and alternative == 'less'):
            edge_cluster = cluster_mass
            if observed:
                edge_cluster_list.append((edge_cluster, inds_list))

    if observed:
        return edge_cluster_list
    else:
        return edge_cluster


def coords_list_to_mask(coords_list, mask_shape):

    mask = np.empty(mask_shape, dtype=bool)
    mask[:, :] = True
    mask[tuple(zip(*coords_list))] = False

    return mask

def vector_angle(v1, v2):

    unit_vector_1 = v1 / np.linalg.norm(v1)
    unit_vector_2 = v2 / np.linalg.norm(v2)
    dot_product = np.dot(unit_vector_1, unit_vector_2)
    angle = np.arccos(dot_product)
    return degrees(angle)
