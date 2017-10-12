import os
import itertools
import pandas as pd
import StringIO
from Bio import Phylo
from pypopgen.modules import haplotools as hap
from pypopgen.modules import vcfpandas as vp
from pypopgen.modules import tree as treelib
import numpy as np
eu = os.path.expanduser
jn = os.path.join

def hap_df_to_tree(hap_df):
    dm = hap.pairwise_diff_mat(hap_df)
    tree = treelib.dm_to_tree(dm)
    tree_str = StringIO.StringIO()
    Phylo.write(tree,tree_str,format='newick')
    return tree_str.getvalue()

def write_genome_trees_hap_vcf(fn,chrom, ana_name, 
                                samples, sample_to_group, snp_window, chunksize=None,
                                path='~'):
    
    def get_df(chunk0,chunk1):
        chunk0.columns = pd.MultiIndex.from_arrays([chunk0.columns,[0]*len(chunk0.columns)])
        chunk1.columns = pd.MultiIndex.from_arrays([chunk1.columns,[1]*len(chunk1.columns)])
        return pd.concat([chunk0,chunk1],axis=1).sortlevel(axis=1)
    

    header, _ = vp.parse_vcf_header(fn)
    vcf_samples = header[9:] 
    
    samples0 = [s[0] for s in samples if s[1] == 0]
    samples1 = [s[0] for s in samples if s[1] == 1]
    
    t0 = vp.get_vcf_df(fn, chrom, header=header, usecols=['CHROM','POS'] + samples0,
                            converters=vp.converters.first_haplotype(samples0), chunksize=chunksize)
    t1 = vp.get_vcf_df(fn, chrom, header=header, usecols=['CHROM','POS'] + samples1,
                            converters=vp.converters.second_haplotype(samples1), chunksize=chunksize)
    
    hap_df = pd.concat(get_df(chunk0,chunk1) for chunk0,chunk1 in itertools.izip(t0,t1))
    try:
        samples_no_data = np.setdiff1d(samples,list(hap_df.columns.values))
    except ValueError,e:
        print samples
        print list(hap_df.columns.values)
        print hap_df
        raise e
    assert not len(samples_no_data), "No data for {}".format(samples_no_data)
    
    focal_hap_df = hap_df.loc[:,samples]
    group_af = focal_hap_df.groupby(sample_to_group, axis=1).mean()
    
    focal_hap_df.columns = [c[0] + '_h' + str(c[1]) for c in focal_hap_df.columns]
    
    #only keep taxon informative
    segregating = ((group_af>0)&(group_af<1)).any(axis=1)
    taxon_informative = (group_af>0).sum(axis=1)>1
    focal_hap_df_ti = focal_hap_df.loc[segregating&taxon_informative,:]
    trees = focal_hap_df_ti.groupby(np.arange(len(focal_hap_df_ti)) // snp_window).apply(hap_df_to_tree)
    locations = focal_hap_df_ti.groupby(np.arange(len(focal_hap_df_ti)) // snp_window)\
                        .apply(lambda df: np.mean(df.index.droplevel(0).values))
    locations.to_csv(jn(eu(path),"tree_locations_{}_{}_snp{}.txt".format(ana_name,
                                                                                   chrom,snp_window)),index=False)
    with open(jn(eu(path),"trees_{}_{}_snp{}.newick".format(ana_name,
                                                                                   chrom,snp_window)),'w') as f:
        for tree in trees:
            f.write(tree)

def get_topos(path, chrom, run_name, snp_window, ngroups=3):
    ntopo_dic = {3:3,4:15}
    assert ngroups in ntopo_dic.keys()
    f = open(jn(path, "weights_{}_{}_snp{}.tsv".format(run_name,chrom,snp_window)))
    topo_strs = [] 
    for i in range(ntopo_dic[ngroups]):
        topo_strs.append(f.readline().strip().split()[1][:-1])
    topos = pd.read_csv(f,sep='\t')
    locs = pd.read_csv(jn(path,"tree_locations_{}_{}_snp{}.txt".format(run_name,chrom,snp_window)),
                       sep='\t',header=None,index_col=0)
    topos.index = locs.index
#   topos.columns = ['introgression','genome','ILS']
    topos.columns = topo_strs
    return topos, topo_strs


