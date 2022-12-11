#!usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Calculate distances between coarse-grained structures
"""

import os
import glob
from pathlib import Path
from itertools import combinations, combinations_with_replacement, product
from multiprocessing import Pool

import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist

from ._PropertyParams import GranthamSim, SM_THREAD_NORM_Sim, PMBEC_Sim
from ._default import DEF_shape_param, DEF_subtype_anchors, DEF_anchor_dir

class CG_SD_Calculator():

    def __init__(self, SimilarityMatrix="Grantham") -> None:
        
        # # fixed parameters
        self.sigma = 0.5
        self.k = 1

        if SimilarityMatrix == "Grantham":
            self.SimMtx = GranthamSim

        elif SimilarityMatrix == "SM_THREAD_NORM":
            self.SimMtx = SM_THREAD_NORM_Sim

        elif SimilarityMatrix == "PMBEC":
            self.SimMtx = PMBEC_Sim

        else:
            raise ValueError(f"Similarity Matrix not recognized: {SimilarityMatrix}")
            
    def ParamExtract(self, DATFile):

        """Extract structure information from CG structure csv file"""

        DAT = pd.read_csv(DATFile)

        resi = DAT['Residue'].values.reshape((-1,1))
        
        # convert back HIS variants
        # resi = resi.reshape(-1).tolist()
        # resi = np.array(resi)

        coord = DAT[['X', 'Y', 'Z']].values
        weight = DAT['Weight'].values.reshape((-1,1))

        accepted = np.flatnonzero(weight) # remove unaccepted residues with weight 0

        return coord[accepted], resi[accepted].reshape(-1).tolist(), weight[accepted]

    def ResiPairSim(self, ResiA, ResiB):
        
        """Residue physicochemical similarity"""
        
        ResiPairComb = [(x, y) for x in ResiA for y in ResiB]
        # print(ResiPairComb)
        ResiPairSim = np.array([self.SimMtx[i][j] for i,j in ResiPairComb])

        return ResiPairSim.reshape((len(ResiA), len(ResiB)))

    def CloudSimilarity(self, CoordA, ResiA, CoordB, ResiB, WeightA, WeightB):
        
        ResiPairSim_score = self.ResiPairSim(ResiA, ResiB)

        SimScore = np.sum( np.reciprocal(np.cosh(self.sigma*cdist(CoordA, CoordB, "euclidean")))**self.k * ResiPairSim_score * np.sqrt(np.outer(WeightA, WeightB)) )

        return SimScore

    def CalcSim(self, DATpair):
        
        """Similarity score between two point clouds"""

        # print(f"Simi: {comb}")
        
        CoordA, ResiA, WeightA = self.ParamExtract(DATpair[0])
        CoordB, ResiB, WeightB = self.ParamExtract(DATpair[1])
        
        return (DATpair, self.CloudSimilarity(CoordA, ResiA, CoordB, ResiB, WeightA, WeightB))
    
    def SaveDist(self, OutCSV):

        """Save distance matrix to csv file"""
        
        self.DistMat.to_csv(OutCSV)
        
        return
    
    def CalcDistMtx(self, DATDir, ListFile=None):
        """
        Pairwise structure distance matrix between PDB files in a directory
        ======================================
        Input:
            DATDir: Directory of coarse-grained HLA structures
            ListFile (optional): List file for selecting alleles, see "../Dataset_split" directory
        """
        DATList = [InDAT for InDAT in glob.glob(f"{DATDir}/*.csv")]
        if ListFile:
            with open(ListFile, "r") as ALF:
                RequestedList = [f"{DATDir}/{line.strip()}" for line in ALF]
                if set(RequestedList)-set(DATList):
                    raise ValueError(f"File in list is not found: {list(set(RequestedList)-set(DATList))}")
                else:
                    DATList = list(set(RequestedList)&set(DATList))
        lite_DATList = [Path(DAT).stem for DAT in DATList] # remove directory name and suffix
        AlleleComb_wi = combinations_with_replacement(DATList, 2)
        AlleleComb_wo = combinations(DATList, 2)

        self.DistMat = pd.DataFrame(np.zeros((len(lite_DATList), len(lite_DATList))), index=lite_DATList, columns=lite_DATList)

        SimilarityMat = {}

        pool = Pool(os.cpu_count())
        result = pool.map(self.CalcSim, AlleleComb_wi)
        
        pool.close()
        pool.join()
        
        for i in result:
            SimilarityMat[i[0]] = i[1]
        
        for comb in AlleleComb_wo:
            # print(f"Dist: {comb}")
            distance = np.sqrt(SimilarityMat[(comb[0], comb[0])] + SimilarityMat[(comb[1], comb[1])] - 2 * SimilarityMat[comb])
            self.DistMat.loc[Path(comb[1]).stem, Path(comb[0]).stem] = distance

        return

    def SelfSim(self, FilePath):
        
        """Similarity Score of the molecule with it self"""

        Coord, Resi, Weight = self.ParamExtract(FilePath)
        
        return (FilePath, self.CloudSimilarity(Coord, Resi, Coord, Resi, Weight, Weight))

    def CalcAnchorDist(self, InDir, anchor_Dir, anchor_dict:dict, ListFile=None):

        """
        Distance of target structures to specified anchor structures
        =======================================
        Input:
            InDir: Directory containing target structures
            anchor_dict: Dictionary of anchor alleles and represented supertype / sub-supertype
            ListFile (optional): Specify target files
            anchor_Dir : Directory containing anchor structures
        """

        TList = [InDAT for InDAT in glob.glob(f"{InDir}/*.csv")] # list of target alleles path
        
        if ListFile:
            with open(ListFile, "r") as ALF:
                RequestedList = [f"{InDir}/{line.strip()}" for line in ALF]
                if set(RequestedList)-set(TList):
                    raise ValueError(f"File in list is not found: {list(set(RequestedList)-set(TList))}")
                else:
                    TList = list(set(RequestedList))

        TList_name = [Path(DAT).stem for DAT in TList] # list of target alleles name

        AList_name = list(set(anchor_dict.keys()))
        AList = [f"{anchor_Dir}/{f}.csv" for f in AList_name]
        if set(AList) - set([DAT for DAT in glob.glob(f"{anchor_Dir}/*.csv")]):
            raise ValueError(f"Anchor allele file(s) not found: {set(AList) - set([DAT for DAT in glob.glob(f'{anchor_Dir}/*.csv')])}")

        self.DistMat = pd.DataFrame(np.zeros((len(TList_name), len(AList_name))), index=TList_name, columns=AList_name)

        AlleleComb = list(product(TList, AList))

        SimilarityMat = {}
        SelfSimMat = {}

        pool = Pool(os.cpu_count())
        inter_result = pool.map(self.CalcSim, AlleleComb)
        intra_result = pool.map(self.SelfSim, AList+TList)
        
        pool.close()
        pool.join()
        
        for i in inter_result:
            SimilarityMat[i[0]] = i[1]

        for i in intra_result:
            SelfSimMat[i[0]] = i[1]

        for comb in AlleleComb:
            # print(f"Dist: {comb}")
            distance = np.sqrt(SelfSimMat[comb[0]] + SelfSimMat[comb[1]] - 2 * SimilarityMat[comb])
            self.DistMat.loc[Path(comb[0]).stem, Path(comb[1]).stem] = distance

        return

# ====================================================================================
def CGCalcMat(CGDATDir, SimMtx, sigma, k, AlleleListFile=None, DistMat_output=None, Standardize:bool=False):
    """
    Calculte pairwise structure distance matrix using coarse-grained distance metric
    ====================================
    Input:
        CGDATDir: Directory of coarse-grained HLA structures
        SimMtx (optional): Similarity matrix, choose from ["Grantham", "SM_THREAD_NORM", "PMBEC"]
        AlleleListFile (optional): List file for selecting alleles, see "../Dataset_split" directory
        sigma, k (optional): Shape parameters
        DistMat_output (optional): File name of distance matrix
        Standardize (optional): If true, standardize the distance matrix to [0-1]

    Output:
        Distance_matrix
    """
    metric = CG_SD_Calculator(SimilarityMatrix=SimMtx)

    metric.sigma = sigma
    metric.k = k

    metric.CalcDistMtx(CGDATDir, AlleleListFile)

    if Standardize:
        metric.DistMat = (metric.DistMat - metric.DistMat.min().min()) / (metric.DistMat.max().max() - metric.DistMat.min().min())

    if DistMat_output:
        metric.SaveDist(DistMat_output)

    return metric.DistMat

def CGAnchorMat(CGDATDir, SimMtx, sigma, k, CGAnchorDir=DEF_anchor_dir, anchor_dict:dict=DEF_subtype_anchors, AlleleListFile=None, DistMat_output=None):
    """
    Calculate distances of query alleles to predifined anchor alleles
    =======================================
    Input:
        CGDATDir: directory of CGDAT files of query alleles
    """
    metric = CG_SD_Calculator(SimilarityMatrix=SimMtx)

    metric.sigma = sigma
    metric.k = k
    
    metric.CalcAnchorDist(CGDATDir, CGAnchorDir, anchor_dict, AlleleListFile)

    if DistMat_output:
        metric.SaveDist(DistMat_output)

    return metric.DistMat