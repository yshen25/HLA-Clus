#!usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Functions for tuning parameters
"""
import pandas as pd
import numpy as np

from .cluster import hierarchical_cluster

def SSE(Mat:pd.DataFrame, clusters:list, square=True):
    # SSE divided by number of elements in each cluster
    if not square:
        Mat = Mat.add(Mat.T, fill_value=0)
    
    sum_SSE = 0

    for clst in clusters:
        
        clst_mat = Mat.loc[clst, clst]
        centroid = clst_mat.sum().idxmin()
        sum_SSE += clst_mat[centroid].pow(2).sum()
    
    return sum_SSE

def Silhouette(Mat:pd.DataFrame, groups:list, square=False):
    if not square:
        Mat = Mat.add(Mat.T, fill_value=0)
    score = []

    if len(groups) == 1:
        return 0

    for i in range(len(groups)):

        if len(groups[i]) == 1:
            # if only one element in a group, the silhouette score is 0 (arbitrary)
            score.append(0)
            continue

        out_groups = groups[0:i] + groups[i+1:]

        for allele in groups[i]:
            
            # distance within groups
            ai = Mat.loc[groups[i],allele].sum() / (len(groups[i]) - 1)
            # distance to neighbor cluster
            bi = np.min([Mat.loc[out_group,allele].sum() / len(out_group) for out_group in out_groups])

            if ai <= bi:
                score.append(1-ai/bi)

            else:
                score.append(bi/ai-1)

        # print(score)

    return np.mean(score)

# parameter tuning
def Tuning_N(ClusterMat, Nmin, Nmax, RefMat=None, ClusterSilhouette=False, RefSilhouette=False) -> tuple:
    """
    Determining optimized number of clusters
    Arguments:
        ClusterMat: Distance matrix that the clustering is based on, lower triangular form
        RefMat: Distance matrix that the performance accessment is based on, lower triangular form
        Nmin: starting number of clusters
        Nmax: ending number of clusters
    
    Output:
        (StructSSE, BASSE, StructSilhouette, BASilhouette)
    """
    ClusterSSE = []
    RefSSE = []
    Silhouette_C = []
    Silhouette_R = []

    # dist_list = []

    for i in range(Nmin, Nmax+1):

        # initialize optional parameters
        # BA_err = 'NA'
        # SilhouetteScore = 'NA'
        # BASilhouetteScore = 'NA'

        cluster, _, _ = hierarchical_cluster(ClusterMat, square=True, N=i, L='complete', threshold=None)
        #complete average single
        groups = [group[1].index.tolist() for group in cluster.groupby(cluster)]
        # print(groups)
        
        ClusterSSE.append(SSE(ClusterMat, groups))

        if RefMat is not None:
            RefSSE.append(SSE(RefMat, groups))

        if ClusterSilhouette:
            Silhouette_C.append(Silhouette(ClusterMat, groups))

        if RefSilhouette:
            Silhouette_R.append(Silhouette(RefMat, groups))

    return ClusterSSE, RefSSE, Silhouette_C, Silhouette_R

def Tuning_shape_param():
    # TODO: implement this!
    return