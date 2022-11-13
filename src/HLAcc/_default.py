#!usr/bin/env python3
# -*- coding: utf-8 -*-
"""
default parameters for import
"""
import os
from pathlib import Path

# DEF_ref_pdb = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'dat', '1i4f_Crown.pdb') # template pdb for trimming and aligning
# DEF_anchor_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'dat', 'Anchor_Alleles')
DEF_ref_pdb = str(Path(__file__).resolve().parent.parent.parent.joinpath('dat/1i4f_Crown.pdb'))
DEF_anchor_dir = str(Path(__file__).resolve().parent.parent.parent.joinpath('dat/Anchor_Alleles'))

DEF_WeightDict = {"A":{63:9.9,67:7.6,116:5.8,9:5.4,97:4.7,152:4.6,167:3.8,156:3.7,74:3.6,70:3.2,80:3.0,171:2.9,45:2.8,77:2.7,76:2.4,114:2.3,99:2.1,95:1.5,158:1.4,24:1.2,7:1.0}}

DEF_shape_param = {"SimMtx":"Grantham", "sigma":0.3, "k":1}

DEF_supertype_anchors = {
    "A01_01":"A01-A03-A66",
    "A02_01":"A02",
    "A02_03":"A02",
    "A02_06":"A02",
    "A02_07":"A02",
    "A03_01":"A01-A03-A66",
    "A11_01":"A01-A03-A66",
    "A24_02":"A24",
    "A30_01":"A01-A03-A66",
    "A68_01":"A02",
    "B07_02":"B07-B35",
    "B08_01":"B08-B18-B39",
    "B14_02":"B14",
    "B15_01":"B15-B40",
    "B18_01":"B08-B18-B39",
    "B27_05":"B27",
    "B35_01":"B07-B35",
    "B39_01":"B08-B18-B39",
    "B40_01":"B15-B40",
    "B40_02":"B15-B40",
    "B42_01":"B07-B35",
    "B44_02":"B44",
    "B44_03":"B44",
    "B46_01":"C01-C02",
    "B51_01":"B51-B58",
    "B57_01":"B51-B58",
    "B58_01":"B51-B58",
    "C04_01":"C01-C02",
    "C05_01":"C01-C02",
    "C06_02":"C01-C02",
    "C08_02":"C01-C02",
    "A26_01":"A01-A03-A66",
    "C07_01":"C07"
}

DEF_subtype_anchors = {
    "A01_01":"A01",
    "A02_01":"A02",
    "A02_03":"A02",
    "A02_06":"A02",
    "A02_07":"A02",
    "A03_01":"A03",
    "A11_01":"A03",
    "A24_02":"A24",
    "A30_01":"A03",
    "A68_01":"A02",
    "B07_02":"B07",
    "B08_01":"B08",
    "B14_02":"B14",
    "B15_01":"B15",
    "B18_01":"B18",
    "B27_05":"B27",
    "B35_01":"B35",
    "B39_01":"B39",
    "B40_01":"B40",
    "B40_02":"B15",
    "B42_01":"B07",
    "B44_02":"B44",
    "B44_03":"B44",
    "B46_01":"C02",
    "B51_01":"B51",
    "B57_01":"B58",
    "B58_01":"B58",
    "C04_01":"C01",
    "C05_01":"C01",
    "C06_02":"C02",
    "C08_02":"C01",
    "A26_01":"A66",
    "C07_01":"C07"
}