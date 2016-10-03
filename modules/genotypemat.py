import gc
import numpy as np
import pandas as pd


def pairwise_diff_numpy(gen_arr):
    """Squared pairwise distances between all 
    columns of 0,1,2 genotype array arr.
    This matrix based function is at least 10 
    times faster than iterating over columns.
    """
    gen_arr = gen_arr.astype(np.float64)-1
    #compare heterozygous with hom alt
    mat1 = np.where(gen_arr==0,-1,gen_arr)
    mat1[np.isnan(gen_arr)]=0
    #and hom ref
    mat2 = np.where(gen_arr==0,1,gen_arr)
    mat2[np.isnan(gen_arr)]=0
    #account for heterozygous comparisons
    mat3 = np.where(gen_arr==0,1,0)
    #don't count nan comparisons
    n = np.dot((~np.isnan(gen_arr)*1.).T,~np.isnan(gen_arr)*1.)
    del gen_arr
    B = (np.dot(mat1.T,mat1)+np.dot(mat2.T,mat2))/2.
    del mat1
    del mat2
    gc.collect()
    het = np.dot(mat3.T,mat3)
    return ((n- B)+het)/2. + np.diag(np.diag(((n- B)+het)/2.)) #self comparisons should not be divided by 2

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


