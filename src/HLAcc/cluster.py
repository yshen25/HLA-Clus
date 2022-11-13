#!usr/bin/env python3
"""
High-level tools for HLA clustering
"""
from sklearn.cluster import DBSCAN, AgglomerativeClustering
import pandas as pd
import numpy as np
from scipy.optimize import linear_sum_assignment

from .visual import plot_dendrogram

# ======= Clustering =======
def hierarchical_cluster(Mat, N=None, L='complete', threshold=None, outtree=None, plot_dendro=False, figsize=(10,5), color_threshold=None, labelsize=12):
    """
    Hierarchical clustering based on square pairwise distance matrix
    ================================
    Input:
        Mat: pairwise distance matrix
        N (optional): number of clusters, one (and only one) of N and threshold must be provided
        square (optional): if input matrix is square matrix (True) or lower triangular matrix (False)
        L (optional): linkage, default = 'complete'
        threshold (optional): distance threshold for a cluster, one (and only one) of N and threshold must be provided
        outtree (optional): file name for output newick tree file
        plot_dendro (optional): if draw dendrogram
        color_threshold (optional): distance threshold for coloring branches as a cluster in dendrogram

    Output:
        (clustering_result, tree_label_order, cluster_centroids)
    """
    if N is None and threshold is None:
        raise ValueError("Both N and threshold not given. At least one need to be specified")

    if not np.allclose(Mat, Mat.T):
        Mat = Mat.add(Mat.T, fill_value=0)

    model = AgglomerativeClustering(n_clusters=N, affinity='precomputed', linkage=L, distance_threshold=threshold, compute_distances=True, compute_full_tree=True).fit(Mat)
    result = pd.Series(model.labels_, index=Mat.index)

    order = None
    
    if plot_dendro:
        order = plot_dendrogram(model, None, Mat.index, figsize, color_threshold, outtree, labelsize)

    centers = []
    for i in result.groupby(by=result):
        group = i[1].index.to_numpy()
        centers.append(Mat.loc[group,group].sum(axis=0).idxmin())

    return result, order, pd.Series(centers, index=range(len(centers)))

def NearestNeighbor_cluster(Mat:pd.DataFrame, anchor_dict:dict) -> pd.DataFrame :
    """
    Nearest neighbor clustering based on distances to anchor alleles
    ================================
    Input:
        Mat: distance matrix to specified anchor structures
        anchor_dict: dictionary specifying anchor alleles and represented supertypes

    Output:
        DataFrame of target allele, nearest anchor allele, and clustering result
    """
    NN_cluster = Mat.T.idxmin().to_frame(name="Nearest_anchor")
    NN_cluster["Cluster"] = NN_cluster["Nearest_anchor"].map(anchor_dict)

    return NN_cluster

def DBSCAN_cluster(Mat:pd.DataFrame, epsilon:float, MinSample=5):
    """
    Cluster alleles using DBSCAN
    In development
    """
    if not np.allclose(Mat, Mat.T):
        Mat = Mat.add(Mat.T, fill_value=0)
    clustering = DBSCAN(eps=epsilon, min_samples=MinSample, metric="precomputed", n_jobs=-1).fit(Mat)
    labels = clustering.labels_
    result = pd.Series(labels, index=Mat.index)

    return result

# ========== stability ===========
def _jaccard(member_index1:list, member_index2:list) -> float:
    # calculate jaccard index between two groups
    s1 = set(member_index1)
    s2 = set(member_index2)
    return len(s1.intersection(s2)) / len(s1.union(s2))

def max_jaccard(clustering1:np.array, clustering2:np.array):
    """
    since cluster name might change, max jaccard is used to indicate matching cluster between
    reference clustering (full sample) and bootstrap clustering
    """
    groups1 = np.unique(clustering1) # name of groups
    groups2 = np.unique(clustering2)
    jar_matrix = np.zeros((len(groups1),len(groups2))) # jaccard index matrix, storing jaccard index of all-to-all groups

    for i in range(len(groups1)):
        for j in range(len(groups2)):
            members1 = np.where(clustering1 == i)[0] # index of group members
            members2 = np.where(clustering2 == j)[0]
            jar_matrix[i,j] = _jaccard(members1, members2)

    row_idx, col_idx = linear_sum_assignment(jar_matrix, maximize=True) # find the solution with maximum sum
    # print(jar_matrix, row_ind, col_ind)
    # print(groups1[row_ind])
    return row_idx, jar_matrix[row_idx, col_idx]

def _bootstrap_sampling(sample_size:int, nb):
    """
    return index of samples
    """
    samples_idx = []
    idx = [i for i in range(sample_size)]
    for _ in range(nb):
        x = np.random.choice(idx, sample_size, replace=True).tolist()
        samples_idx.append(x)

    return samples_idx

def cluster_stability(dist_mat:pd.DataFrame, ref_clustering:pd.Series, NB:int=100, average:bool=True)->np.array:
    """
    returns average/step-wise jaccard index (similarity) of each cluster
    """
    # distance matrix is upper-triangle, change to square form
    dist_mat = dist_mat.add(dist_mat.T, fill_value=0)
    dist_mat = dist_mat.to_numpy()

    ref_groups = np.unique(ref_clustering)
    N_groups = len(ref_groups)
    N_samples = dist_mat.shape[0]

    jac_matrix = np.empty((NB, N_groups))
    jac_matrix[:] = np.nan
    bt_index = _bootstrap_sampling(N_samples, NB)
    for i in range(NB):
        index = bt_index[i]
        dist_mat_bt = dist_mat[:,index][index,:]
        dist_mat_bt = pd.DataFrame(dist_mat_bt)
        ref_group = ref_clustering[index]

        bt_group, _, _ = hierarchical_cluster(dist_mat_bt, N=N_groups)

        group_idx, group_jac = max_jaccard(ref_group, bt_group)

        jac_matrix[i, group_idx] = group_jac

    if average:
        return pd.Series(np.nanmean(jac_matrix, axis=0))

    else:
        return jac_matrix