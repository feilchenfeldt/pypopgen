import os, gc
import numpy as np
import pandas as pd

def pairwise_diff_mat(df):
    """
    Calculate pairwise difference data frame.
    For Genotype 0,1,2 data frame.
    Uses numpy matrix multiplication.
    """
    cols = df.columns
    diff = pairwise_diff_numpy(df.values)
    gc.collect()
    return pd.DataFrame(diff,index=cols,columns=cols)

def pairwise_diff_numpy(gen_arr):
    """
    ATTENTION: This does not handle nans
    for missing data implement a mask 
    Setting missing to zero beforehand would be wrong.
    """
    n = np.dot(gen_arr.T,np.logical_not(gen_arr)) + \
        np.dot(np.logical_not(gen_arr.T),gen_arr)
    return n

def pw_diff_to_individual(hap_pairwise_diff):
    ind_pairwise_diff = th.sum(axis=1,level=0).sum(axis=0, level=0)
    denum =  np.ones(ind_pairwise_diff.shape) * 4.
    np.fill_diagonal(denum, 2.)
    return ind_pairwise_diff / denum
