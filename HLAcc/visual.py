#!usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Functions for visulizing results and evaluating performance
"""
from scipy.cluster.hierarchy import dendrogram, to_tree
import seaborn as sn
from scipy.stats import linregress

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def dist_heatmap(Mat, order:list=None, size=(10,10), label=False, line=False, labelsize=8, **cbar_kw):
    """
    Visualize distance matrix as heatmap
    ======================================
    Input:
        Mat: distance matrix
        order: list of allele names. Order of heatmap is re-arranged according to the list. If line == True, order need to be nested list
        size: size of heatmap
        label: show/hide allele names
        line: show/hide solid line between clusters. The clusters are difined by nested list in order, for example: [[A,B],[C,D]]
        labelsize: font size of label
    """

    if not np.allclose(Mat, Mat.T):
        Mat = Mat.add(Mat.T, fill_value=0)

    if order:
        if type(order[0]) == list:
            flat_order = [item for sublist in order for item in sublist]
        else:
            flat_order = order
        Mat = Mat[flat_order] # re-arrange row order
        Mat = Mat.reindex(flat_order) # re-arrange column order

    # standardize to [0,1]
    Mat = (Mat - Mat.min().min()) / (Mat.max().max() - Mat.min().min())
    
    plt.figure(figsize=size)
    if label:
        ticks = True
    else:
        ticks = False

    g = sn.heatmap(Mat, square=True, xticklabels=ticks, yticklabels=ticks, cbar_kws=cbar_kw)
    g.axes.tick_params(axis='both', labelsize=labelsize, pad=10)

    # draw lines between clusters
    if line:
        split = np.cumsum([len(sublist) for sublist in order])
        for line in split[:-1]:
            plt.axhline(y=line, color='k', linestyle='-')
            plt.axvline(x=line, color='k', linestyle='-')
    plt.show()
    return

def anchor_heatmap(Mat, order:list=None, size=(10,5), label=False, line=False, labelsize=8, **cbar_kw):
    if order:
        if type(order[0]) == list:
            flat_order = [item for sublist in order for item in sublist]
        else:
            flat_order = order
        Mat = Mat[flat_order] # re-arrange row order
        Mat = Mat.reindex(flat_order) # re-arrange column order

    # standardize to [0,1]
    Mat = (Mat - Mat.min().min()) / (Mat.max().max() - Mat.min().min())
    # vertical to horizontal
    Mat = Mat.T
    
    plt.figure(figsize=size)
    if label:
        ticks = True
    else:
        ticks = False

    g = sn.heatmap(Mat, xticklabels=ticks, yticklabels=ticks, cbar_kws=cbar_kw)
    g.axes.tick_params(axis='both', labelsize=labelsize, pad=10)

    # draw lines between clusters
    if line:
        split = np.cumsum([len(sublist) for sublist in order])
        for line in split[:-1]:
            plt.axhline(y=line, color='k', linestyle='-')
            plt.axvline(x=line, color='k', linestyle='-')
    plt.show()
    return

def _getNewick(node, newick, parentdist, leaf_names):
    """
    convert dendrogram to newick format
    """
    if node.is_leaf():
        return "%s:%.2f%s" % (leaf_names[node.id], parentdist - node.dist, newick)
    else:
        if len(newick) > 0:
            newick = "):%.2f%s" % (parentdist - node.dist, newick)
        else:
            newick = ");"
        newick = _getNewick(node.get_left(), newick, node.dist, leaf_names)
        newick = _getNewick(node.get_right(), ",%s" % (newick), node.dist, leaf_names)
        newick = "(%s" % (newick)
        return newick

def plot_dendrogram(model, truncate, labels, figsize, color_threshold, outtree=None, labelsize=12):
    """
    Plot the dendrogram
    ==================================
    Input:
        model
        truncate
        labels
        color_threshold
    """
    # create the counts of samples under each node
    counts = np.zeros(model.children_.shape[0])
    n_samples = len(model.labels_)
    for i, merge in enumerate(model.children_):
        current_count = 0
        for child_idx in merge:
            if child_idx < n_samples:
                current_count += 1  # leaf node
            else:
                current_count += counts[child_idx - n_samples]
        counts[i] = current_count

    linkage_matrix = np.column_stack([model.children_, model.distances_, counts]).astype(float)

    # Plot the corresponding dendrogram
    plt.figure(figsize=figsize)

    dendro = dendrogram(linkage_matrix, truncate_mode=truncate, labels=labels, leaf_font_size=labelsize, get_leaves=True, color_threshold=color_threshold)

    plt.show()

    if outtree:
        tree = to_tree(linkage_matrix)
        OutFile = _getNewick(tree, "", tree.dist, labels)
        with open(outtree, "w") as fh:
            fh.write(OutFile)

    return dendro['leaves']

def correlation(ArrayA, ArrayB, show_plot=True, xlabel="SD", ylabel="PD"):
    """
    Correlation plot between two array
    =============================
    Input:
        ArrayA, ArrayB: 1-D array with same shape, A as x and B as y
        show_plot: if true, draw correlation plot
    
    Output:
        (slope, intercept, rvalue)
    """

    # standardize to 0-1
    xx = (ArrayA - np.min(ArrayA)) / (np.max(ArrayA) - np.min(ArrayA))
    yy = (ArrayB - np.min(ArrayB)) / (np.max(ArrayB) - np.min(ArrayB)) 

    slope, intercept, rvalue, _, _ = linregress(xx, yy)

    if intercept >= 0:
        equation = f"Y = {round(slope, 2)}X+{round(intercept, 2)}"
    else:
        equation = f"Y = {round(slope, 2)}X{round(intercept, 2)}"
    # intercept = 0
    if show_plot:
        plt.figure(figsize=(8,8))
        plt.xlim(0,1)
        plt.ylim(0,1)
        plt.scatter(xx, yy, )

        x_vals = np.array([0, 1])
        y_vals = intercept + slope * x_vals
        # y_vals = slope * x_vals
        plt.plot(x_vals, y_vals, '--', linewidth=3, c='k')
        plt.xticks(fontsize=26)
        plt.yticks(fontsize=26)
        plt.xlabel(xlabel, fontsize=28)
        plt.ylabel(ylabel, fontsize=28)

        plt.text(0.02,0.95, equation, fontsize=30)
        plt.text(0.02,0.89, f"R = {round(rvalue, 2)}", fontsize=30)
        plt.show()
    
    return (slope, intercept, rvalue)

def elbow_plot(ClusterSSE, Silhouette_C, RefSSE=None, Silhouette_R=None, Nmin=1, Nmax=12):
    """
    Draw elbow plot (Sum-of-Squared-Error (SSE) versus number of clusters (N)) to determine number of clusters
    """

    xx = range(Nmin, Nmax+1)

    lines = []
    labels = []

    fig, ax1 = plt.subplots(figsize=(6,10))
    ax2 = ax1.twinx()

    line1, = ax1.plot(xx, ClusterSSE, c='b', marker='^', mfc='None', mec='b', ms='8', mew=3, alpha=0.6, label="SSE")
    line2, = ax2.plot(xx, Silhouette_C, c='gold', marker='v', mfc='None', mec='g', ms='8', mew=3, alpha=0.6, label="Silhouette")
    lines.append(line1)
    labels.append("SSE")
    lines.append(line2)
    labels.append("Silhouette")

    if RefSSE is None and Silhouette_R is None:
        pass
    elif None in (RefSSE, Silhouette_R):
        raise ValueError("Both RefSSE and Silhouette_R must be provided or leave None")
    else:
        line3, = ax1.plot(xx, RefSSE, c='b', marker='^', mfc='None', mec='b', ms='8', mew=3, alpha=0.6, label="Ref SSE")
        line4, = ax2.plot(xx, Silhouette_R, c='g', marker='v', mfc='None', mec='g', ms='8', mew=3, alpha=0.6, label="Ref Silhouette")

        lines.append(line3)
        labels.append("Ref SSE")
        lines.append(line4)
        labels.append("Ref Silhouette")

    ax1.set_xlabel('Number of clusters', fontsize=20)
    ax1.set_xticks(range(1,Nmax+1,2))
    ax1.tick_params(axis='x', labelsize=16)
    
    ax1.set_ylabel('SSE', color='tab:blue', fontsize=20)
    ax1.tick_params(axis='y', labelcolor='tab:blue', labelsize=16)

    ax2.set_ylabel('Silhouette', color='tab:green', fontsize=20)
    ax2.tick_params(axis='y', labelcolor='tab:green', labelsize=16)

    ax1.legend(lines, labels, prop={"size":16})
    ax1.grid(linestyle='--')

    # fig.legend()
    plt.show()

    return