import os, csv
import itertools
import pandas as pd
import StringIO
from Bio import Phylo

from pypopgen.modules import vcfpandas as vp
from pypopgen.modules import treetools
from pypopgen.modules import tensorfstats

from pypopgen.modules import vcfpandas as vp


#REMOVE THESE DEPENDENCIES IN THE FUTURE
from pypopgen.modules import haplotools as hap
from pypopgen.modules import tree as treelib

import numpy as np
eu = os.path.expanduser
jn = os.path.join


def write_genome_trees(vcf_fn,chrom, ana_name,
                                samples, snp_window, chunksize=None,
                                path='~'):
    hap_df = vp.get_haplotype_df(vcf_fn, chrom,
                   #     start=1e7, end=1e7+100000,
                       samples=samples)
    segregating = (hap_df>0).any(axis=1)&(hap_df<1).any(axis=1)
    missing = hap_df.isnull().any(axis=1)
    hap_df1 = hap_df[segregating&(~missing)]
    window_grps = hap_df1.groupby(np.arange(len(hap_df1)) // snp_window)
    names = [i[0] for i in hap_df1.columns.values[::2]]
    newicks = window_grps.apply(lambda df: 
                treetools.dm_to_tree(
                    tensorfstats.calc.divergence(
                            df.groupby(axis=1,level=0)),
                                    names=names).write(format=5))
    locations = window_grps.apply(lambda df: np.mean(df.index.droplevel(0).values))
    locations.to_csv(jn(eu(path),"tree_locations_{}_{}_snp{}.txt".format(ana_name,
                                                                chrom,snp_window)),index=False)
    pd.DataFrame(newicks).to_csv(jn(eu(path),"trees_{}_{}_snp{}.newick".format(ana_name,
                                                                chrom,snp_window)),
               index=False, header=False, quoting=csv.QUOTE_NONE, sep='\t')

def get_topos(chrom,ana_name,snp_window,ntopos=3, path='~'):
    ntopo_dic = {3:3, 4:15}
    assert ntopos in ntopo_dic.keys(), "{} not in {}".format(ntopos, ntopo_dic.keys())
    f = open(jn(path, "weights_{}_{}_snp{}.tsv".format(ana_name,chrom,snp_window)))
    topo_strs = [] 
    for i in range(ntopo_dic[ntopos]):
        topo_strs.append(f.readline().strip().split()[1][:-1])
    topos = pd.read_csv(f,sep='\t')
    locs = pd.read_csv(jn(path,"tree_locations_{}_{}_snp{}.txt".format(ana_name,
                                                                                  chrom,snp_window)),
                       sep='\t',header=None,
                       index_col=0)
    
    topos.index = locs.index
#   topos.columns = ['introgression','genome','ILS']
    topos.columns = topo_strs
    return topos, topo_strs

def get_topos_proc(chrom,ana_name,snp_window, ntopos=3 ,path='~'):
    topos, tstrl = get_topos(chrom,ana_name,snp_window, ntopos=ntopos, path=path)
    topos = topos.div(topos.sum(axis=1),axis=0)
    return topos

def get_topos_all_chrom(chromosomes,ana_name,snp_window, ntopos=3, path='~'):
    twist_df = []
    for chrom in chromosomes:
        topos_df = get_topos_proc(chrom,ana_name=ana_name, snp_window=snp_window, ntopos=ntopos, path=path)
        topos_df.index = pd.MultiIndex.from_product([chrom,topos_df.index])
        topos_df.index.names = ['chrom','pos']
        twist_df.append(topos_df)
    twist_df = pd.concat(twist_df)
    return twist_df


#___________________________________________________________
#DEPRECIATED": THE FOLLOWING FUNCTIONS SHOULD NOT BE USED
#
#def hap_df_to_tree(hap_df):
#    dm = hap.pairwise_diff_mat(hap_df)
#    tree = treelib.dm_to_tree(dm)
#    tree_str = StringIO.StringIO()
#    Phylo.write(tree,tree_str,format='newick')
#    return tree_str.getvalue()
#
#def write_genome_trees_hap_vcf(fn,chrom, ana_name, 
#                                samples, sample_to_group, snp_window, chunksize=None,
#                                path='~'):
#    
#    def get_df(chunk0,chunk1):
#        chunk0.columns = pd.MultiIndex.from_arrays([chunk0.columns,[0]*len(chunk0.columns)])
#        chunk1.columns = pd.MultiIndex.from_arrays([chunk1.columns,[1]*len(chunk1.columns)])
#        return pd.concat([chunk0,chunk1],axis=1).sortlevel(axis=1)
#    
#
#    header, _ = vp.parse_vcf_header(fn)
#    vcf_samples = header[9:] 
#    
#    samples0 = [s[0] for s in samples if s[1] == 0]
#    samples1 = [s[0] for s in samples if s[1] == 1]
#    
#    t0 = vp.get_vcf_df(fn, chrom, header=header, usecols=['CHROM','POS'] + samples0,
#                            converters=vp.converters.first_haplotype(samples0), chunksize=chunksize)
#    t1 = vp.get_vcf_df(fn, chrom, header=header, usecols=['CHROM','POS'] + samples1,
#                            converters=vp.converters.second_haplotype(samples1), chunksize=chunksize)
#    
#    hap_df = pd.concat(get_df(chunk0,chunk1) for chunk0,chunk1 in itertools.izip(t0,t1))
#    try:
#        samples_no_data = np.setdiff1d(samples,list(hap_df.columns.values))
#    except ValueError,e:
#        print samples
#        print list(hap_df.columns.values)
#        print hap_df
#        raise e
#    assert not len(samples_no_data), "No data for {}".format(samples_no_data)
#    
#    focal_hap_df = hap_df.loc[:,samples]
#    group_af = focal_hap_df.groupby(sample_to_group, axis=1).mean()
#    
#    focal_hap_df.columns = [c[0] + '_h' + str(c[1]) for c in focal_hap_df.columns]
#    
#    #only keep taxon informative
#    segregating = ((group_af>0)&(group_af<1)).any(axis=1)
#    taxon_informative = (group_af>0).sum(axis=1)>1
#    focal_hap_df_ti = focal_hap_df.loc[segregating&taxon_informative,:]
#    trees = focal_hap_df_ti.groupby(np.arange(len(focal_hap_df_ti)) // snp_window).apply(hap_df_to_tree)
#    locations = focal_hap_df_ti.groupby(np.arange(len(focal_hap_df_ti)) // snp_window)\
#                        .apply(lambda df: np.mean(df.index.droplevel(0).values))
#    locations.to_csv(jn(eu(path),"tree_locations_{}_{}_snp{}.txt".format(ana_name,
#                                                                                   chrom,snp_window)),index=False)
#    with open(jn(eu(path),"trees_{}_{}_snp{}.newick".format(ana_name,
#                                                                                   chrom,snp_window)),'w') as f:
#        for tree in trees:
#            f.write(tree)



