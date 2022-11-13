#!usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd

def square2triangle(InMat):
    keep = np.invert(np.triu(np.ones(InMat.shape)).astype('bool'))
    return InMat.mask(keep, other=0)

def triangle2square(InMat):
    return InMat.add(InMat.T, fill_value=0)

def crop_mtx(Mtx:pd.DataFrame, order:list, flatten:bool=False):
    
    if type(order[0]) == list:
        flat_order = [item for sublist in order for item in sublist]
    else:
        flat_order = order
    
    Mtx = Mtx.loc[flat_order, flat_order] # re-arrange row order

    if flatten:
        # since matrix is symmetric, each pairwise distance exist 2 copies
        # first extract lower triangle (no diagonal), then extract values to make sure only extract unique values
        keep = np.invert(np.triu(np.ones(Mtx.shape)).astype('bool')).flatten()
        return Mtx.to_numpy().flatten()[keep]

    return Mtx