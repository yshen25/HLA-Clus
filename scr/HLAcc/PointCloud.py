#!usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Extract atom and residue coordinates from protein structures in pdb format
The atomic or coarse-grained residue coordinates are stored in csv files
"""
import os
import glob
from string import digits
from itertools import combinations
from pathlib import Path
from tempfile import NamedTemporaryFile
from typing import Union
from itertools import groupby
from operator import itemgetter
import shutil

import numpy as np
# from scipy.stats import hypsecant
import pandas as pd
from pymol import cmd
from Bio.Align import PairwiseAligner
from Bio.PDB.DSSP import dssp_dict_from_pdb_file

from biopandas.pdb import PandasPdb
import matplotlib.cm as cm
# import open3d as o3d

from .PropertyParams import AtomicMass, PartialCharge, AtomicHydrophobicity

def _Seqinfo_from_struct(Struct):
    # Struct = PandasPdb().read_pdb(PDBfile)
    Seq_df = Struct.amino3to1()

    info = {} # chain: (start, stop, sequence)
    for group in Struct.df['ATOM'].groupby(['chain_id']):
        # group[0]: chain id
        sequence = ''.join(list(Seq_df[Seq_df['chain_id'] == group[0]]['residue_name']))
        info[group[0]] = [group[1].iloc[0]['residue_number'], group[1].iloc[-1]['residue_number'], sequence]

    return info

def _trim_renumber_struct(Struct, Qinfo, Tinfo):
    trimmed = []
    for chain in Qinfo:
        Qstart = Qinfo[chain][0]
        Qend = Qinfo[chain][1]
        trimmed_chain = Struct.df['ATOM'].query(f"chain_id == '{chain}' and {Qstart} <= residue_number <= {Qend}").copy()

        Tstart = Tinfo[chain][0]
        shift = Tstart - Qstart # how many positions are shifted
        trimmed_chain.loc[:,'residue_number'] += shift

        trimmed.append(trimmed_chain)
    return pd.concat(trimmed)

def PDB_trim(InDir, TemplatePDB, OutDir):
    # template struct object and chain info
    Tstruct = PandasPdb().read_pdb(TemplatePDB)
    Tinfo = _Seqinfo_from_struct(Tstruct)

    aligner = PairwiseAligner()
    aligner.gap_score = -10 # no gap wanted
    aligner.target_left_open_gap_score = -5

    for InPDB in glob.glob(f"{InDir}/*.pdb"):
        # query struct and chain info
        Qstruct = PandasPdb().read_pdb(InPDB)
        Qinfo = _Seqinfo_from_struct(Qstruct)

        # remove unwanted chains
        for key in Qinfo:
            if key not in Tinfo:
                Qinfo.pop(key, None)

        for chain in Tinfo:
            Tseq = Tinfo[chain][2]
            Qseq = Qinfo[chain][2]
            alignments = aligner.align(Qseq, Tseq)
            # start query position that matched with template
            qstart = alignments[0].aligned[0][0][0] # 0 based and relative, the first resi is always 0
            # qend = alignments[0].aligned[0][-1][-1]

            # Qend
            # Qinfo[chain][1] = Qinfo[chain][0] + qend # 1 based and indexed
            Qinfo[chain][1] = qstart + Tinfo[chain][1] - Tinfo[chain][0] + 1 # fixed length
            # Qstart
            Qinfo[chain][0] = Qinfo[chain][0] + qstart
            
        Qstruct.df['ATOM'] = _trim_renumber_struct(Qstruct, Qinfo, Tinfo)# trimed as marked in query info, renumbered as the template
        
        Qstruct.to_pdb(f"{OutDir}/{Path(InPDB).name}", records=['ATOM'])
    
    return

def BB_RMSD(Qpdb, Rpdb):

    cmd.load(Qpdb, "query")
    cmd.load(Rpdb, "ref")

    rmsd = cmd.align(f"query////CA", f"ref////CA", cycles=0, transform=0)[0] # backbone RMSD

    cmd.delete("query")
    cmd.delete("ref")

    return rmsd

def _centroid_structure(InPDBDir):
    """
    find the structure that has least RMSD with other structures
    used as alignment reference
    """
    # list all pdb files in input dir
    InFiles = glob.glob(f"{InPDBDir}/*.pdb")
    
    # initialize RMSD matrix
    RMSD_Mat = pd.DataFrame(np.zeros((len(InFiles), len(InFiles)), dtype=float), index=InFiles, columns=InFiles)

    # calculate pairwise RMSD and fill in the matrix
    comb = combinations(InFiles, 2)
    for pair in comb:
        RMSD_Mat.loc[pair[1],pair[0]] = BB_RMSD(pair[1], pair[0])

    RMSD_Mat = RMSD_Mat.add(RMSD_Mat.T, fill_value=0) # lower triangle to square
    # print(RMSD_Mat)
    RMSD_sum = RMSD_Mat.sum(axis=1)
    return RMSD_sum.idxmin()

def _alphaNbeta(InPDB, dssp_path=None):
    """
    Extract alpha helix and beta sheets aa index of input PDB (1-based)
    Return string of index range
    ====================
    InPDB: Input PDB file
    dssp_path: Path to mkdssp excutable
    --------------------
    output
    anchor_index: Indexes of residues that forms alpha helix and beta strands
    anchor_range: Ranges of indexes
    """
    if dssp_path:
        if not os.path.exists(dssp_path):
            raise ValueError("DSSP not found!")
    else:
        dssp_path = shutil.which("mkdssp")
        if not dssp_path:
            raise ValueError("DSSP not found!")
    
    dssp = dssp_dict_from_pdb_file(InPDB, DSSP=dssp_path)
    secondary_structure = [dssp[0][i][1] for i in dssp[0].keys()]
    aa_index = [dssp[0][i][5] for i in dssp[0].keys()] # 1-based

    anchor_index = np.array(aa_index)[np.isin(secondary_structure, ["E", "H", "I", "G"])]

    anchor_range = []
    for _, g in groupby(enumerate(anchor_index),lambda x:x[0]-x[1]):
        group = list(map(itemgetter(1),g))
        anchor_range.append(f"{group[0]}-{group[-1]}")

    anchor_range = "+".join(anchor_range)

    return anchor_index, anchor_range

def PDB_align(InDir, OutDir, refPDB=None, SSAlign=False, AddH_only=False):
    """
    Superimpose input PDBs to the reference PDB
    if reference PDB not specified and use geomatric center alignment, the centroid structure is used
    Hydrogens are added for Glycines
    if AddH_only, no alignment is performed
    """
    def change_resname(df, resnum, changeto):
        df.loc[df["residue_number"] == resnum, "residue_name"] = changeto
        return df

    def change_atomname(df, index, changeto):
        df.at[index, "atom_name"] = changeto
        return df

    if not os.path.exists(OutDir):
        print(f"Create directory for aligned PDBs: {OutDir}")
        os.mkdir(OutDir)
    
    if refPDB is None:
        if AddH_only:
            refPDB = glob.glob(f"{InDir}/*.pdb")[0]

        else:
            # Find the centroid structure as reference
            print("Reference structure is not specified\n\tSearching for centroid structure...")
            refPDB = _centroid_structure(InDir)
            print(f"Centroid: {refPDB}")

    if SSAlign:
        _, RAchRange = _alphaNbeta(refPDB)

    cmd.load(refPDB, "template")
    tmp = NamedTemporaryFile().name

    for InPDB in glob.glob(f"{InDir}/*.pdb"):
        # =============== PyMOL: align and add H ================
        # print("align:", InPDB)
        cmd.load(InPDB, "target")
        cmd.h_add(selection="(resn 'GLY' and name 'CA')") # add hydrogen atoms for Glycine, especifically for crystal structures
        cmd.alter("(resn 'GLY' and name 'H01')", "name='1HA'")
        cmd.alter("(resn 'GLY' and name 'H02')", "name='2HA'")
        cmd.alter("(resn 'GLY' and name 'H03')", "name='2HA'")

        if AddH_only:
            pass

        elif SSAlign: # align and superimpose based on secondary structures
            _, TAchRange = _alphaNbeta(InPDB)
            cmd.align(f"target///{TAchRange}/CA", f"template///{RAchRange}/CA")
        
        else:
            cmd.align(f"target////CA", f"template////CA") # align and superimpose

        cmd.save(tmp, "target",0,'pdb')
        cmd.delete("target")

        # ================ biopandas: correct residue and atom name ==============
        Qstruct = PandasPdb().read_pdb(tmp)
        new_df = Qstruct.df['ATOM']
        for residue in Qstruct.df['ATOM'].groupby('residue_number'):

            if residue[1].iloc[0]["residue_name"] == "HIS":

                if "HE2" in residue[1]['atom_name'].values:
                    if "HD1" in residue[1]['atom_name'].values:
                        change_resname(new_df, residue[0], "HIP") # H on both epsilon and delta N
                    else:
                        change_resname(new_df, residue[0], "HIE") # H on epsilon N only

                else:
                    change_resname(new_df, residue[0], "HID") # H on delta N only. By default HIS is HID

            elif residue[1].iloc[0]["residue_name"] == "MSE":
                change_resname(new_df, residue[0], "MET") # no partial charge value for MSE

        # now residue is the last residue
        idx = np.where(residue[1]['atom_name'].values == 'OXT')[0]
        if idx:
            change_atomname(new_df, residue[1].iloc[idx].index[0], "O") # Rosetta relax will automaticly change last Oxigen from O to OXT

        Qstruct.df['ATOM'] = new_df
        Qstruct.to_pdb(f"{OutDir}/{Path(InPDB).name}", records=["ATOM"])
        
    cmd.delete("template")
        
    return

def PDB_to_AtomCloud(InDir, OutDir):
    """
    Assign parameters like partial charge to each atom, and store dataframe into csv file
    """
    # with open('ffcharge.pkl', 'rb') as inf:
    #     ffcharge = pickle.load(inf)

    ffcharge = PartialCharge
    ffhydro = AtomicHydrophobicity

    # with open('pep_surf.pkl', 'rb') as inf:
    #     pep = pickle.load(inf)

    if not os.path.exists(OutDir):
        os.makedirs(OutDir)

    for InPDB in glob.glob(f"{InDir}/*.pdb"):

        pro = PandasPdb().read_pdb(InPDB)
        ATOMS = pro.df["ATOM"]
        # atom_coord = ATOMS[['x_coord', 'y_coord', 'z_coord']]
        # ligand_coord = pro.df['HETATM'][['x_coord', 'y_coord', 'z_coord']]

        # atom to ligand distance is the minimum of atom to ligand-atom distances
        # dist2ligand = cdist(atom_coord, ligand_coord, 'euclidean')
        # depth = np.min(dist2ligand, axis=1)

        # depth and acceptance are set to 1 initially
        depth = accpt = np.ones((ATOMS.shape[0],1))
        
        # atomic partial charge and hydrophobicity
        basic_info = ATOMS[["residue_name", "chain_id", "residue_number", "atom_name", "atom_number", "x_coord", "y_coord", "z_coord"]].to_numpy()
        charge_hydro = []
        for row in basic_info:
            atom = row[3].lstrip(digits)
            resi = row[0]

            charge = ffcharge[resi][atom]

            if resi in ["HIE", "HIP", "HID"]: # change back to normal HIS
                row[0] = resi = "HIS"

            hydro = np.nan if atom.startswith("H") else ffhydro[resi][atom] # hydrophobicity doesn't distinguish HIS types
            
            
            charge_hydro.append([charge, hydro])

        OutList = np.hstack([basic_info, charge_hydro, depth, accpt])

        OutDF = pd.DataFrame(OutList, columns=["Residue", "Chain", "ResNum", "Atom", "AtomNum", "X", "Y", "Z", "Charge", "Hydrophobicity", "Depth", "Accept"])

        OutDF.to_csv(f"{OutDir}/{Path(InPDB).stem}.csv", index=False)

    return

def CoarseGrain(InDAT, OutCGDAT):
    """
    Coarse graining of single csv file
    center-of-mass of side chain
    """
    # if OutFile is None:
    #     OutFile = DATFile.split(".")[0] + "_CG.csv"

    FullAtom_df = pd.read_csv(InDAT)
    # DAT = DAT[~DAT.Atom.isin(["N", "CA", "C", "O"])] # filter out backbone atoms

    ResGen = FullAtom_df.groupby(by=['Chain', 'ResNum'])

    CG_row = []
    for resi in ResGen:
        x = y = z = 0
        total_mass = 0
        chain = resi[0][0]
        resnum = resi[0][1]
        resatoms_df = resi[1]
        # print(resatoms)
        resname = resatoms_df["Residue"].iloc[0]

        for atom in resatoms_df[["Atom", "X", "Y", "Z"]].to_numpy():
            # print(atom)
            if atom[0] in ["N", "CA", "C", "O"]:
                continue # filter out backbone atoms

            if not atom[0][0].isdigit():
                mass = AtomicMass[atom[0][0]] # real atom name is the first letter of pdb atom name

            else:
                mass = AtomicMass[atom[0][1]]

            total_mass += mass
            x += mass * atom[1]
            y += mass * atom[2]
            z += mass * atom[3]

        if total_mass == 0:
            CG_row.append([chain, resnum, resname, np.nan, np.nan, np.nan, 0])
        else:
            CG_row.append([chain, resnum, resname, x/total_mass, y/total_mass, z/total_mass, 1])

    CG_DAT = pd.DataFrame(CG_row, columns=["Chain", "ResNum", "Residue", "X", "Y", "Z", "Weight"])
    CG_DAT.to_csv(OutCGDAT, index=False)

    return

def PointCloudCG(DATDir, OutDir):
    """
    Coarse graining of atom cloud files
    ==================================
    Input: full atom DAT directory
    Output: coarse-grained DAT file
    """
    for InDAT in glob.glob(f"{DATDir}/*.csv"):
        CoarseGrain(InDAT, f"{OutDir}/{Path(InDAT).name}")

    return

def element_depth():
    return

def depth_threshold():
    return

def min_dist():
    return

def dist_threshold():
    return

def PDBtoPointCloud(InputDir:Union[str, bytes, os.PathLike], AlignDir:Union[str, bytes, os.PathLike], 
    PCDir:Union[str, bytes, os.PathLike], TrimDir:Union[str, bytes, os.PathLike]=None, CGDir:Union[str, bytes, os.PathLike]=None, 
    SSAlign:bool=False, RefPDB:Union[str, bytes, os.PathLike]=None)->None:
    """
    Turn PDB file into atom cloud
    ==============================
    Input: Input directory that contains PDB files
    SSAlign: Do structure alignment based on residues that form secondary structures, alphahelixes and beta strands
    trim: Trim every input PDB files according to Reference PDB, the Reference must be specified
    RefPDB: PDB file used as reference in structure alignment. If not specified, the centroid structure will be calculated and used
    --------------------
    Output: Trimmed and aligned PDB files, CSV pointcloud files
    """
    # ===== Trim =====
    # If trim, the RefPDB must be specified
    if TrimDir:
        print("Start trimming...")
        if not RefPDB:
            raise ValueError("If enable trim, a RefPDB must be specified")
        if not os.path.exists(TrimDir):
            os.mkdir(TrimDir)
            print(f"\tCreate directory for trimmed PDBs: {TrimDir}")
        PDB_trim(InputDir, RefPDB, TrimDir)
        print(f"\t...Trimming done, saved files to {TrimDir}")
    
    # ===== Align =====
    print("Start Aligning...")
    if not os.path.exists(AlignDir):
        print(f"\tCreate directory for aligned PDBs: {AlignDir}")
        os.mkdir(AlignDir)

    if TrimDir:
        PDB_align(TrimDir, AlignDir, refPDB=RefPDB, SSAlign=SSAlign)
    else:
        PDB_align(InputDir, AlignDir, refPDB=RefPDB, SSAlign=SSAlign)
    print(f"\t...Aligning done, saved files to {AlignDir}")

    # ===== To point cloud =====
    print("Start converting to atom cloud...")
    if not os.path.exists(PCDir):
        print(f"\tCreate directory for atom cloud CSVs: {PCDir}")
        os.mkdir(PCDir)
    PDB_to_AtomCloud(AlignDir, PCDir)
    print(f"\t...Convert to atom cloud done, saved files to {PCDir}")

    # ===== Coarse graining =====
    if CGDir:
        print("Start coarse graining...")
        if not os.path.exists(CGDir):
            print(f"\tCreate directory for coarse-grained point cloud CSVs: {CGDir}")
            os.mkdir(CGDir)
        PointCloudCG(PCDir, CGDir)
        print(f"\t...Coarse graining done, saved files to {CGDir}")
    return

class assign_weight:
    def read_CG_DAT(self, InCG_DAT):
        self.filename = InCG_DAT
        df = pd.read_csv(InCG_DAT)
        self.df = df.set_index(['Chain', 'ResNum'])
        return

    def by_resi_number(self, weight_dict):
        reform = [[i,j,weight_dict[i][j]] for i in weight_dict.keys() for j in weight_dict[i].keys()]
        self.weight = pd.DataFrame(reform, columns=["Chain", "ResNum", "Weight"])
        self.weight = self.weight.set_index(['Chain', 'ResNum'])
        
        return
    def update_weight(self):
        # initialize by giving 0 in every position
        self.df['Weight'] = 0
        
        # residue weight as multi-indexed dataframe
        # self.weight_dict_to_df(weight_dict)

        # change weight via update
        self.df.update(self.weight)
        return

    def save_CG_DAT(self):
        self.df.to_csv(self.filename)
        return

# def Show_PointCloud(DATDir):
#     pointcloudlist = []
#     cmap = cm.ScalarMappable(cmap="coolwarm")
#     for InDAT in glob.glob(f"{DATDir}/*.csv"):
#         df = pd.read_csv(InDAT)
#         P1 = o3d.geometry.PointCloud()
#         coord = df[['X','Y','Z']].to_numpy()
#         weight = df['Weight'].to_numpy()

#         color1 = cmap.to_rgba(weight)[:,0:3]
#         P1.points = o3d.utility.Vector3dVector(coord)
#         P1.colors = o3d.utility.Vector3dVector(color1)
#         pointcloudlist.append(P1)

#     o3d.visualization.draw_geometries(pointcloudlist)
#     return

def reweight_by_dict(InDir, WeightDict):
    """
    change residue weight in coarse grained csv file according to provided dict
    """
    re_weight = assign_weight()
    re_weight.by_resi_number(WeightDict)

    for CGDAT_file in glob.glob(f"{InDir}/*.csv"):

        re_weight.read_CG_DAT(CGDAT_file)
        re_weight.update_weight()
        re_weight.save_CG_DAT()
    return