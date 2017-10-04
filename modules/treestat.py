"""
This module provides tools to investigate genetic relatedness
of multiple samples. It contains methods to investigate
how well the genetic relatedness between samples is represented
by a (phylogenetic) tree. In particular, if there was genetic
exchange (gene flow, introgression, admixture) between non-sister
clades (species of populations), then this would lead to non-tree-like
genetic relatedness. Genetic relatedness would be better described by
a graph. However, this package does not attempt to infer this graph,
but it provides methods to infer the which samples or which internal
branches violate the tree model (and thus likely were subject to
introgression.

Supported input files:
Block compressed and indexed Variant Call Format (.vcf.gz) files
that were compressed with bgzip and indexed with Tabix (.tbi file present).
Note: Things might be faster when working with an intermediate file format
      (Haplotype/genotype matrix) we might consider implementing conversion
      to such intermediate files in the future.
"""
__version__ = '0.1'
__author__ = 'Hannes Svardal'

# global imports
import logging
import pathos.multiprocessing as mp
import pandas as pd
import numpy as np

import matplotlib as mpl
from matplotlib import pyplot as plt


# local imports
import haplotools as hap

logger = logging.getLogger()
logging.basicConfig(format='%(levelname)-8s %(asctime)s %(filename) '
                    '%(message)s')
logger.setLevel(logging.WARNING)


# -----------------------------------
# -----Get pairwise differences------
# -----------------------------------

def get_tree_residual_mat(tree, pwd):
    names_tree_order = [c.name for c in tree.get_terminals()]
    tree_dist = pd.DataFrame(index=names_tree_order,columns=names_tree_order)

    for n0 in tree_dist.index:
        for n1 in tree_dist.columns:
            if n0!=n1:
                tree_dist.loc[n0,n1] = tree.distance(n0,n1)
            else: 
                tree_dist.loc[n0,n1] = np.nan
    tree_dist1 = tree_dist.stack(dropna=False)
    tree_dist1.name = 'tree_dist'
    pw_diff = pwd.stack()
    pw_diff.name = 'pw_diff'
    
    tree_dist_vs_pwd = pd.DataFrame(tree_dist1).join(pw_diff,how='inner')
    discrepancy_pct = (tree_dist_vs_pwd['pw_diff']-tree_dist_vs_pwd['tree_dist'])/tree_dist_vs_pwd['pw_diff']*100           
    discrepancy_mat = discrepancy_pct.unstack()
    return discrepancy_mat.ix[names_tree_order,names_tree_order].astype(np.float).fillna(0)


def plot_residuals(residuals, tree, pwd):
    tree_order = [c.name for c in tree.get_terminals()]

    upper_mat = pwd.ix[tree_order,tree_order]
    upper = upper_mat.copy()
    upper.iloc[:,:] = np.tril(upper.values)
    upper_masked = np.ma.array(upper, mask=(upper==0))

    lower_mat = residuals.ix[tree_order,tree_order].astype(float)
    lower = lower_mat.copy()
    lower.iloc[:,:] = np.triu(lower.values,k=1)
    lower_masked = np.ma.array(lower, mask=(lower==0))


    mpl.rcParams["font.size"]=14
    per_country = False#country labels (else subspecies labels) 

    assert (upper_mat.index == lower_mat.index).all()
    assert (upper_mat.columns == lower_mat.columns).all()
    assert (upper_mat.columns == upper_mat.index).all()

    tree_ax_width = 0.3
    cb_ax_width = 0.05
    cb_width = 0.03

    fig = plt.figure(figsize=(12*(1+tree_ax_width)+2,12+1))



    cb_ax_width_x = cb_ax_width * (1-tree_ax_width)
    cb_width_x = cb_width * (1-tree_ax_width)



    dmax = lower_mat.max().max()
    dmin = lower_mat.min().min()
    colors2 = plt.cm.Reds(np.linspace(0, dmax/abs(dmin), int(256*dmax)))
        
    colors1 = plt.cm.Blues_r(np.linspace(0., 1, int(256*abs(dmin))))
    # combine them and build a new colormap
    colors = np.vstack((colors1, colors2))
    mymap = mpl.colors.LinearSegmentedColormap.from_list('my_colormap', colors)



    tick_locations = []
    tick_names = []
    borders = []
    pos = 0 




    # for g in upper_mat.index.droplevel(1).unique():      
    #     d = len(upper_mat.loc[g].loc[:,g])
    #     #for plotting, call C. trimaculatus a sanddweller
    #     if g == 'utaka':
    #         d=d-1
    #     elif g == 'sand':
    #         d=d+1
    #     tick_names.append(g)
    #     borders.append(pos+d)
    #     tick_locations.append(pos+d/2.)
    #     pos += d

    ax = fig.add_axes([tree_ax_width, 0, 1-tree_ax_width-cb_ax_width_x, 1-cb_ax_width])   

    diff_mesh = ax.pcolormesh(lower_masked,cmap=mymap,rasterized=True)

    cbaxes = fig.add_axes([1-cb_width_x,0, cb_width_x , 1-cb_ax_width])   
    diff_bar = plt.colorbar(diff_mesh,cax = cbaxes,
                                    label='(Pairwise difference - tree distance) / pairwise differences (%)')#,use_gridspec=False,orientation='horizontal',anchor=(0.5,1))
    diff_bar.solids.set_rasterized(True) 


    #cbaxes.tickparams()
    #cbaxes.yaxis.tick_left()




     
    cNorm  = mpl.colors.Normalize(vmin=upper_mat.min().min(), vmax=upper_mat.max().max())
    pi_mesh = ax.pcolormesh(upper_masked,cmap='Greens',norm=cNorm,rasterized=True)#,norm=cNorm

    cbaxes_pi = fig.add_axes([tree_ax_width,1-cb_width, 1-tree_ax_width-cb_ax_width, cb_width]) 
    pi_bar = plt.colorbar(pi_mesh,cax = cbaxes_pi,orientation='horizontal',
                           label=r'Pairwise difference ($\mathrm{kb}^{-1}$)')#,use_gridspec=False,orientation='horizontal',anchor=(0.5,1))
    pi_bar.solids.set_rasterized(True) 
    cbaxes_pi.tick_params(labelsize=16) 
    cbaxes_pi.xaxis.labelpad = 10
    cbaxes_pi.xaxis.set_label_position('top') 
    cbaxes_pi.xaxis.tick_top()




    ax.set_xlim([0,len(upper_mat)])
    ax.set_ylim([0,len(upper_mat)])
    #ax.set_xticks(tick_locations)
    #ax.set_xticklabels(tick_names)#rotation='vertical',ha='center'
    ax.set_xticks(np.arange(len(upper_mat))+0.5)
    labels = upper_mat.index.values


    ax.set_xticklabels(labels,
                                fontsize=9,rotation='vertical',ha='center',fontstyle='normal')
     #ax.set_xticklabels([])



    ax.set_yticks(np.arange(len(upper_mat))+0.5)
    ax.set_yticklabels(labels,fontsize=9,fontstyle='normal')
    #ax.set_yticklabels(tick_names[1:],rotation='vertical')
    #ax.set_yticklabels(dm_group_rad.index.droplevel(0).values,fontsize=8)
    #ax.set_yticks(np.arange(len(dm_group_rad)))
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    ax.tick_params(axis='y', size=0)
    ax.tick_params(axis='x', size=0)
    for b in borders:   
        l = ax.axvline(x=b, ymin=0,ymax=1,linewidth=1., color='w')
        l = ax.axhline(y=b, xmin=0, xmax=1,linewidth=1.,color='w')
                     
    xls, yls = (ax.get_xticklabels(),ax.get_yticklabels())
                     #for i,(group,_) in enumerate(upper_mat.index.values):
                     #    xls[i].set_color(group_colors[group])
                     #    yls[i].set_color(group_colors[group])
                         
                             
                             #tree_ax = fig.add_axes([0., 0, tree_ax_width, 1-cb_ax_width])




                             #dstat.draw_tree(tree_no, axes=tree_ax, do_show=False,
                             #               label_func=lambda c:'')#get_short_scientific(c.name)

    #tree_ax.set_ylim((ax.get_ylim()[0]-0.5,ax.get_ylim()[1]-0.5))
    #tree_ax.axis('off')
    return fig

def get_pairwise_differences(fn, samples=None, ncpus='auto', chunksize=10000):
    """
    Get matrix of pairwise genetic differences.

    fn ... input filename (.vcf or .vcf.gz)
    samples ... Sample ids to use as given in the
                vcf header. If 'None', all samples are used
    ncpus ... Number of processes to spawn. If
              'auto', number of available cpus is used.
    """
    hap.get_pairwise_diff(fn, samples=samples, chunksize=chunksize,
                          map_fun=pool.map_async, get_result_fun=lambda r: r.get())


def get_fstat_snpwindow_chrom(chrom, callset, path, populations,
                              quadruples, controlsamples_h3=1,
                              controlsamples_h2=0, jackknife_window=2000):
    """
    Jackknife window is in number of informative SNPs.
    """
    gen_df = get_gen_df(callset, chrom, path)

    # fstats, jackknife_window_fstats = get_fstat_snpwindow(gen_df, quadruples,
    #                                                       populations,
    #                                                       controlsamples_h3=\
    #                                                       controlsamples_h3,
    #                                                       controlsamples_h2=\
    #                                                       controlsamples_h2,
    #                                                       jackknife_window=\
    #                                                       jackknife_window)
    # return fstats, jackknife_window_fstats
    return gen_df
