import os, gc
import itertools
import numpy as np
import pandas as pd
#local import
import vcfpandas as vp



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
    ind_pairwise_diff = hap_pairwise_diff.sum(axis=1,level=0).sum(axis=0, level=0)
    denum =  np.ones(ind_pairwise_diff.shape) * 4.
    np.fill_diagonal(denum, 2.)
    return ind_pairwise_diff / denum

def get_pairwise_diff_region(fn,samples, chrom,start,end,read_chunksize=50000):
    """
    Read haplotypes from specific region tabixed vcf.gz file
    and calculate pairwise differences.
    """
    def get_pdw(chunk0,chunk1):
        chunk0.columns = pd.MultiIndex.from_arrays([chunk0.columns,[0]*len(chunk0.columns)])
        chunk1.columns = pd.MultiIndex.from_arrays([chunk1.columns,[1]*len(chunk1.columns)])
        hap_chunk =  pd.concat([chunk0,chunk1],axis=1).sortlevel(axis=1)
        dm = pairwise_diff_mat(hap_chunk)
        return dm
    
    samples = [str(s) for s in samples]

    header, _ = vp.parse_vcf_header(fn)
    
    t0 = vp.get_vcf_df(fn, chrom, start=start,end=end,header=header, usecols=['CHROM','POS'] + samples,
                        converters=vp.converters.first_haplotype(samples), chunksize=read_chunksize)
    t1 = vp.get_vcf_df(fn, chrom, start=start,end=end, header=header, usecols=['CHROM','POS'] + samples,
                        converters=vp.converters.second_haplotype(samples), chunksize=read_chunksize)

    dm = reduce(lambda a,b:a+b,(get_pdw(chunk0,chunk1) for chunk0,chunk1 in itertools.izip(t0,t1)))

    dm_ind = pw_diff_to_individual(dm)
    
    return dm_ind
