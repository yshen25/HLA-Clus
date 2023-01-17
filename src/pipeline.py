

import os
from pathlib import Path
from tempfile import TemporaryDirectory

from ._CG_metric import CGAnchorMat, CGCalcMat
from ._default import (DEF_ref_pdb, DEF_shape_param, DEF_anchor_dir, DEF_subtype_anchors,
                       DEF_WeightDict)
from .cluster import (NearestNeighbor_cluster, cluster_stability,
                      hierarchical_cluster)
from .coarsegrain import PDB2PointCloud, reweight_by_dict
from .tools import triangle2square
from .visual import anchor_heatmap, dist_heatmap


def _dir_exist(dirname):
    raise ValueError(f"{dirname} already exist, cannot make directory. Please rename or delete it.")

def Processing_pipeline(PDBDir, CGDir, WeightDict=DEF_WeightDict, clean=False):
    """
    Fully automatic structure processing that turn PDB files into CG point clouds
    All setting are default
    """
    if clean:
        trim_dir = TemporaryDirectory().name
        align_dir = TemporaryDirectory().name
        pointcloud_dir = TemporaryDirectory().name
    else:
        wkDir = str(Path(PDBDir).parent.absolute())

        trim_dir = wkDir + "/TRIM"
        if os.path.exists(trim_dir):
            _dir_exist(trim_dir)

        align_dir = wkDir + "/ALIGN"
        if os.path.exists(align_dir):
            _dir_exist(align_dir)

        pointcloud_dir = wkDir + "/FAcloud"
        if os.path.exists(pointcloud_dir):
            _dir_exist(pointcloud_dir)

    PDB2PointCloud(InputDir=PDBDir, TrimDir=trim_dir, AlignDir=align_dir, PCDir=pointcloud_dir, CGDir=CGDir, RefPDB=DEF_ref_pdb, SSAlign=True)
    reweight_by_dict(CGDir, WeightDict)

    return

# ========== clustering pipeline ===========
def HC_pipeline(CGDATDir, NumClusters, SimMtx:str=None, sigma:float=None, k:float=None, AlleleListFile=None, DistMatName:str=None, DendroName:str=None, HMSize=(5,5), HMLabelsize=10, DendroSize=(10,5), DendroLabelsize=12):
    if DistMatName is None:
        DistMatName = os.getcwd() + "/Distance_Matrix.csv"

    if SimMtx is None:
        SimMtx = DEF_shape_param['SimMtx']
    if sigma is None:
        sigma=DEF_shape_param['sigma']
    if k is None:
        k=DEF_shape_param['k']

    print(f"Applied shape parameters: Similarity Matrix: {SimMtx}, sigma: {sigma}, k: {k}")

    Mat = CGCalcMat(CGDATDir, SimMtx, sigma, k, AlleleListFile=AlleleListFile, DistMat_output=DistMatName)
    Mat = triangle2square(Mat)
    dist_heatmap(Mat, label=True, size=HMSize, labelsize=HMLabelsize)
    supertype_cluster, _, supertype_centroids = hierarchical_cluster(Mat, N=NumClusters, L='complete', plot_dendro=True, figsize=DendroSize, outtree=DendroName, labelsize=DendroLabelsize)
    
    ClusterResult = supertype_cluster.to_frame(name='cluster').sort_values(by='cluster')
    ClusterResult.to_csv(os.getcwd() + "/Cluster_Result.csv")
    print("=========Cluster Result==========")
    print(ClusterResult)
    
    stability = cluster_stability(Mat, supertype_cluster)
    
    # cluster centroids and boot strap stability measured by mean Jaccard index
    StabResult = supertype_centroids.to_frame(name='centroids').join(stability.to_frame(name='stability'))
    StabResult.to_csv(os.getcwd() + "/Cluster_Stability.csv")
    print("=========CLuster Centroids and Stability==========")
    print(StabResult)

    print(f"Distance matrix saved as: {DistMatName}\nClustering result saved as: {os.getcwd() + '/Cluster_Result.csv'}\nCluster centroids and stability saved in: {os.getcwd() + '/Cluster_Stability.csv'}")
    return

def NN_pipeline(CGDATDir, SimMtx:str=None, sigma:float=None, k:float=None, anchor_dir=DEF_anchor_dir, anchor_dict=DEF_subtype_anchors, AlleleListFile=None, DistMatName:str=None, HMSize=(10,5), HMLabelsize=10):
    
    if DistMatName is None:
        DistMatName = os.getcwd() + "/NNDistance_Matrix.csv"

    if SimMtx is None:
        SimMtx = DEF_shape_param['SimMtx']
    if sigma is None:
        sigma=DEF_shape_param['sigma']
    if k is None:
        k=DEF_shape_param['k']

    print(f"Applied shape parameters: Similarity Matrix: {SimMtx}, sigma: {sigma}, k: {k}")

    Mat = CGAnchorMat(CGDATDir, SimMtx, sigma, k, anchor_dir, anchor_dict, AlleleListFile=AlleleListFile, DistMat_output=DistMatName)
    anchor_heatmap(Mat, label=True, size=HMSize, labelsize=HMLabelsize)

    ClusterResult = NearestNeighbor_cluster(Mat, anchor_dict)
    ClusterResult.to_csv(os.getcwd() + "/NNCluster_Result.csv")
    print(ClusterResult)
    
    return