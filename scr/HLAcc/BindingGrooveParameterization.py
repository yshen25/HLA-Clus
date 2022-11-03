#!usr/bin/env python3
# -*- coding: utf-8 -*-

"""
extract coordinates of resi forming the binding groove of HLA moleules
input: pdb file
output: csv file
"""

import os
import re
# import pickle
import shutil
from string import digits
from itertools import chain, groupby
from operator import itemgetter

import numpy as np

from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import PPBuilder
from Bio.Align import PairwiseAligner
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB.ResidueDepth import get_surface
from Bio.PDB.DSSP import dssp_dict_from_pdb_file
# from Bio.PDB import extract

from scipy.spatial import Delaunay
from scipy.spatial.distance import cdist

from pymol import cmd

import pandas as pd

class ChainSelector:
    """
    Adapted from Bio.PDB.Dice module
    Only accepts residues with right chainid, between start and end.
    Remove waters and ligands. Only use model 0 by default.
    Hydrogens are kept
    """

    def __init__(self, chain_id, start, end, model_id=0):
        """Initialize the class."""
        self.chain_id = chain_id
        self.start = start
        self.end = end
        self.model_id = model_id

    def accept_model(self, model):
        """Verify if model match the model identifier."""
        # model - only keep model 0
        if model.get_id() == self.model_id:
            return 1
        return 0

    def accept_chain(self, chain):
        """Verify if chain match chain identifier."""
        if chain.get_id() in self.chain_id:
            return 1
        return 0

    def accept_residue(self, residue):
        """Verify if a residue sequence is between the start and end sequence."""
        # residue - between start and end
        hetatm_flag, loc, icode = residue.get_id()
        chain_ = residue.parent.get_id()
        if hetatm_flag != " ":
            # skip HETATMS
            return 0
        if icode != " ":
            print(f"WARNING: Icode {icode} at position {loc}")
        if self.start[self.chain_id.index(chain_)] <= loc <= self.end[self.chain_id.index(chain_)]:
            return 1
        return 0

    def accept_atom(self, atom):
        """Modified to accept all atoms including hydrogens"""
        
        # _hydrogen = re.compile("[123 ]*H.*")

        # if _hydrogen.match(atom.get_id()): # remove hydrogens
        #     return 0
        if "H" in atom.get_id(): # new way to remove hydrogens
            return 0

        if atom.altloc not in [" ", "A"]: # remove altloc atoms
            return 0

        return 1

def extract(structure, chain_id, start, end, filename):
    """Write out selected portion of structure to <filename (pdb)>."""
    sel = ChainSelector(chain_id, start, end)
    io = PDBIO()
    io.set_structure(structure)
    io.save(filename, sel)

    return

def PDB_renumber(struct, start:list):
    """
    renumber pdb files, the first residue in pdb file is indexed as start number
    each chain is given a start number in alphabat order
    """
    for i, chain in enumerate(struct[0]):
    #loop 1: renumber residues to negative number to avoid errors
        residue_id = -1
        for residue in chain.get_residues():
            residue.id=(' ',residue_id,' ')
            residue_id -= 1
        #loop 2
        residue_id = start[i]
        for residue in chain.get_residues():
            #print(chain.get_id(), residue_id)
            residue.id=(' ',residue_id,' ')
            residue_id += 1
    return struct

def PDB_trim(InDir, TemplatePDB, OutDir, OutCSV, chain="A", length=[179], template_start_id=[2]):
    """
    PDB structure trim to have same length with tamplate
    """
    # length=[179] # HLA1
    # start_id=[2] # HLA1
    # length = [80, 85] # HLA2
    # template_start_id = [4,7] # HLA2
    record = []

    PepBuilder = PPBuilder()
    parser = PDBParser(PERMISSIVE=True, QUIET=True)

    TStruct = parser.get_structure("template", TemplatePDB)
    TSeqs = PepBuilder.build_peptides(TStruct)
    
    aligner = PairwiseAligner()
    aligner.gap_score = -10 # no gap wanted

    for InPDB in os.listdir(InDir):
        if InPDB.endswith(".pdb"):
            print("trim:", InPDB)
            InStruct = parser.get_structure("target", f"{InDir}/{InPDB}")
            Seqs = PepBuilder.build_peptides(InStruct[0])
            qends = []
            qstarts = []
            for i, chain_identifier in enumerate(chain):
                InSeq = Seqs[i].get_sequence()
                TSeq = TSeqs[i].get_sequence()
            
                starting_loc = InStruct[0][chain_identifier].child_list[0]._id[1] # index of the first residue
                alignments = aligner.align(InSeq, TSeq)

            ## === archive === stand-alone trim
            # qstart = alignments[0].aligned[0][0][0] + starting_loc # alignment is 0-based, starting loc is 1-based
            # qend = qstart + 178 # fixed length
            #qend = alignments[0].aligned[0][-1][-1] # aligned portion
            # use 177 to remove last amino acid of relaxed models
            
            ## === use with PDB_renumber ===
                # qstart = starting_loc - alignments[0].aligned[0][0][0] + template_start_id[i] -1
                qstart = starting_loc + alignments[0].aligned[1][0][0] - alignments[0].aligned[0][0][0] + template_start_id[i] -1
                # starting loc is calibarated, for the 1st residue in template is not loc 1
                qend = length[i]

                qstarts.append(qstart)
                qends.append(qend)
                
                record.append([InPDB, chain_identifier, qstart, qend, qend-qstart+1])

            # OutPDB = InPDB.split("S")[0].replace("*", "").replace(":", "_") + ".pdb"
            OutPDB = InPDB

            # InStruct = PDB_renumber(InStruct, qstart)
            # extract(InStruct, chain, start_id, qends, f"{OutDir}/{OutPDB}") # HLA1
            InStruct = PDB_renumber(InStruct, qstarts)
            extract(InStruct, chain, [2,2], qends, f"{OutDir}/{OutPDB}")
            #print(f"Trim file saved: {OutDir}/{OutPDB}, {qend-qstart+1}")
            

    df = pd.DataFrame(record, columns=["FileName", "chain", "qstart", "qend", "length"])
    df.to_csv(OutCSV)

    return

def alphaNbeta(InPDB):
    """
    Extract alpha helix and beta sheets aa index of input PDB (1-based)
    Return string of index range
    """
    if os.path.exists("/home/shawn/local/dssp/mkdssp"):
        DSSP_path = "/home/shawn/local/dssp/mkdssp"
    elif os.path.exists("/Users/ys0/local/dssp/mkdssp"):
        DSSP_path = "/Users/ys0/local/dssp/mkdssp"
    else:
        raise ValueError("DSSP not found!")
    
    dssp = dssp_dict_from_pdb_file(InPDB, DSSP=DSSP_path)
    secondary_structure = [dssp[0][i][1] for i in dssp[0].keys()]
    aa_index = [dssp[0][i][5] for i in dssp[0].keys()] # 1-based

    anchor_index = np.array(aa_index)[np.isin(secondary_structure, ["E", "H", "I", "G"])]

    anchor_range = []
    for _, g in groupby(enumerate(anchor_index),lambda x:x[0]-x[1]):
        group = list(map(itemgetter(1),g))
        anchor_range.append(f"{group[0]}-{group[-1]}")

    anchor_range = "+".join(anchor_range)

    return anchor_index, anchor_range

def PDB_align(InDir, refPDB, OutDir):
    """
    Superimpose query PDB to template PDB
    """
    
    cmd.load(refPDB, "template")
    _, RAchRange = alphaNbeta(refPDB)

    for InPDB in os.listdir(InDir):
        if InPDB.endswith(".pdb"):
            print("align:", InPDB)
            _, TAchRange = alphaNbeta(f"{InDir}/{InPDB}")

            cmd.load(f"{InDir}/{InPDB}", "target")
            cmd.h_add(selection="(resn 'GLY' and name 'CA')") # add hydrogen atoms for Glycine, especifically for CG methods of crystal structures
            cmd.alter("(resn 'GLY' and name 'H01')", "name='1HA'")
            cmd.alter("(resn 'GLY' and name 'H02')", "name='2HA'")
            cmd.alter("(resn 'GLY' and name 'H03')", "name='2HA'")
            cmd.align(f"target///{TAchRange}/CA", f"template///{RAchRange}/CA") # align and superimpose based on alpha helix wall and beta sheet plate

            OutPDB = f"{OutDir}/{InPDB.split('.')[0]}.pdb"
            cmd.save(OutPDB, "target")
            # print(f"Align file saved: {OutPDB}")
            cmd.delete("target")
    
    return

def groove_CA_coord(PDBpath, Struct):
    """
    Input pdb and corresponding model parsed by PDBparser
    return array of coordinates of CA that form helix wall and sheet plate
    """
    anchor_index, _ = alphaNbeta(PDBpath)
    OutList = []
    for chain in Struct:
        for residue in chain:
            if residue.id[1] in anchor_index:
                atom = residue["CA"]
                X_coord, Y_coord, Z_coord = atom.coord[0:3]
                OutList.append([X_coord, Y_coord, Z_coord])
    
    return np.array(OutList)

def in_groove(PDBpath, Struct, AtomCoord, PepSurf):
    """
    input pdb, corresponding model, and atom coordinate
    return a 0/1 array of atom inside/outside of binding groove defined by helix wall and sheet plate
    """
    box = Delaunay(groove_CA_coord(PDBpath, Struct))
    InOrOut = box.find_simplex(AtomCoord)>=0
    InOrOut = InOrOut.astype(int).reshape(-1,1)

    NearPep = np.min(cdist(AtomCoord, PepSurf),1) <= 5.0
    NearPep = NearPep.astype(int).reshape(-1,1)

    return InOrOut * NearPep

def resi_depth(Struct, AtomCoord, PepSurf):
    """
    first extact the binding groove surface, then calculate the distance of each atom to binding groove
    """
    # whole protein surface
    surface = get_surface(Struct, MSMS="/home/shawn/local/msms_i86_64Linux2_2.6.1/msms.x86_64Linux2.2.6.1")
    hull = Delaunay(AtomCoord)

    # remove surface vertices that fall out side of convex hull
    in_pocket = hull.find_simplex(surface)>=0
    pocket = surface[in_pocket,:]

    # remove surface vertices that is far from bound peptide
    in_groove = np.min(cdist(pocket, PepSurf),1) <= 1.0 # controls distance cut-off between peptide surface and binding groove surface
    groove = pocket[in_groove,:]

    ResiDepth = np.min(cdist(AtomCoord, groove),1).reshape(-1,1)

    return ResiDepth

def PDB_preprocess(PDBDIr, TemplatePDB, TrimDir, AlignDir, OutCSV, **kwargs):

    if not os.path.exists(TrimDir):
        os.makedirs(TrimDir)

    if not os.path.exists(AlignDir):
        os.makedirs(AlignDir)

    PDB_trim(PDBDIr, TemplatePDB, TrimDir, OutCSV, **kwargs)
    PDB_align(TrimDir, TemplatePDB, AlignDir)

    return

if __name__ == "__main__":

    ## ====models====
    PDB_preprocess("../HLA1_models/CF_relaxed/PDB", "1i4f_Crown.pdb", "../HLA1_models/CF_relaxed/TRIM", "../HLA1_models/CF_relaxed/ALIGN", "HLA1_CF_relaxed_trim.csv")

    pass