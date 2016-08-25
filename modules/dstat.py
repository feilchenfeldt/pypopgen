import os, copy, itertools
import numpy as np
import pandas as pd
import matplotlib as mpl
from matplotlib import pyplot as plt


eu = os.path.expanduser
jn = os.path.join

def get_gen_df(callset, chrom, path):
    df =  pd.read_csv(jn(path,"{}_genotypes_012_{}.tsv".format(callset,chrom)), 
                       sep='\t',index_col=[0,1], na_values='N')
    #dropna to remove sites where the outgoups have no allele
    return df.dropna()



#Functions for tree pruning and search of (h1,h2,h3,o) consistent with the populations

def get_consistent_quadruples(populations, tree, outgroups) :
    consistent_quadruples = []
    
    try:
        pops = populations.keys()
    except AttributeError:
        pops = populations
        
    for tpl in itertools.combinations(pops,4):
        for h1,h2,h3,o in itertools.permutations(tpl):
            try:
                if consistent_with_tree2(tree, h1,h2,h3,o) and \
                    (h1,h2,h3,o) not in consistent_quadruples and o in outgroups:
                    consistent_quadruples.append((h2,h1,h3,o)) 
            except ValueError, e:
                print (h1,h2,h3,o)
                raise e
    return consistent_quadruples

def consistent_with_tree2(in_tree, h1,h2,h3,o):

    tree = copy.deepcopy(in_tree)
    for tc in in_tree.get_terminals():
        if tc.name not in [h1,h2,h3,o]:
            tree.prune(target=tc.name)
    if get_n_nodes(tree,h1,h2) < get_n_nodes(tree,h2,h3) and \
               get_n_nodes(tree,h1,h2) < get_n_nodes(tree,h1,h3) and \
               get_n_nodes(tree,h2,h3) < get_n_nodes(tree,h2,o) and \
               get_n_nodes(tree,h1,h3) < get_n_nodes(tree,h1,o):
        return True
    else:
        return False


def get_n_nodes(tree, leaf1, leaf2):
    return len(tree.trace(leaf1,leaf2))-1

def prune_tree_to_populations(in_tree, population_dict, prune_missing=True):
    tree = copy.deepcopy(in_tree)
    for pop, names in population_dict.iteritems():
        try:
            #TODO: make branch length mean of all subclades!?
            [i for i in tree.find_clades(names[0])][0].name = pop
            
        except IndexError:
            print "could not find clade", names[0]
            print  [i for i in tree.find_clades(names[0])]
        for n in names[1:]:
            tree.prune(target=n)
    if prune_missing:
        for clade in tree.get_terminals():
            if clade.name not in population_dict.keys():
                tree.prune(clade)
    return tree

#------------------------------------------------------------------
#Functions to calculate D-statistics, bootstrap and D over windows
#------------------------------------------------------------------

def get_dstat_snpwindow_chrom(chrom, callset, path, populations, quadruples, jackknife_window=2000, snp_window=None):
    """
    Jackknife window is in number of informative SNPs.
    """
    gen_df = get_gen_df(callset, chrom, path)
    pop_s = pd.Series(
                {n:k for k,v in populations.iteritems() for n in v})
    g = gen_df.groupby(pop_s,axis=1)
    af = g.mean()/2.
    #del g
    #del gen_df
    #gc.collect()
    dstats, jackknife_window_dstats, snp_window_dstats = get_dstat_snpwindow(af, quadruples, 
                                                                     jackknife_window=2000, snp_window=snp_window)
    return dstats, jackknife_window_dstats, snp_window_dstats

def get_dstat_snpwindow(af, quadruples, jackknife_window=2000, snp_window=None):
    #min fraction of snps to report value 
    #(only makes a difference for right-most interval)
    min_observation_fraction=0.75
    dstats = []
    jackknife_window_dstats = []
    snp_window_dstats = []
    for j,(h1,h2,h3,o) in enumerate(quadruples):
        dstat = pd.DataFrame(columns=['num','denom'])
        dstat['num'] = ((af[h1] - af[h2])*(af[h3] - af[o])).dropna()
        dstat['denom']= ((af[h1] + af[h2] - 2 * af[h1] * af[h2])*(af[h3] + af[o] - 2 * af[h3] * af[o])).dropna()
        #only use informative SNPs
        dstat = dstat[dstat['denom'] != 0]
        dstats.append([dstat['num'].sum(),dstat['denom'].sum()])
        jackknife_window_sum = pd.rolling_sum(dstat,jackknife_window,
                   min_periods=int(min_observation_fraction*jackknife_window),
                       center=True).iloc[jackknife_window/2::jackknife_window].dropna()
        jackknife_window_dstats.append(jackknife_window_sum.reset_index(level=1).values.tolist())
        #del jackknife_window_sum
        #del dstat
        #gc.collect
        if snp_window is not None:
            snp_window_sum = pd.rolling_sum(dstat,snp_window,
                   min_periods=min_observation_fraction*snp_window,
                       center=True).iloc[snp_window/2::snp_window].dropna()
            #gc.collect()
            snp_window_dstats.append(snp_window_sum.reset_index(level=1).values.tolist())
    return dstats,jackknife_window_dstats, snp_window_dstats


def dstat_from_array(num_demom_arr):
    return np.sum(num_demom_arr[:,0])/np.sum(num_demom_arr[:,1])


def reduce_dstat_snpwindow(result, chromosomes=None):
    
    def add_chrom(i,tpl_ix):
        chrom = chromosomes[i]
        res = result[i][2][tpl_ix]
        chromvec = np.tile(chrom,len(res))
        chromvec.shape = (chromvec.shape[0],1)
        return np.hstack((chromvec,res))
    
    d = []
    Z = []
    window_d = []
    #iterate over index of h1,h2,h3,o tuples
    for tpl_ix in range(len(result[0][0])):

        nums = [result[i][0][tpl_ix][0] for i in range(len(result))]
        denoms = [result[i][0][tpl_ix][1] for i in range(len(result))]
        d.append(sum(nums)*1./sum(denoms))

        jackknife_arr = np.concatenate([result[i][1][tpl_ix] for i in range(len(result)) if result[i][1][tpl_ix]])[:,1:]

        dstat_from_jackknifes = dstat_from_array(jackknife_arr)

        jackknife_estimates = [dstat_from_array(jackknife_arr[np.arange(jackknife_arr.shape[0])!=i]) \
                                                           for i in range(jackknife_arr.shape[0])]
        print "Number of jackknife estimates for quadruple {}:".format(tpl_ix), len(jackknife_estimates)
        
        Z.append(dstat_from_jackknifes/(np.std(jackknife_estimates,ddof=1)*np.sqrt(len(jackknife_estimates))))

        if result[0][2]:
            window_arr = np.concatenate([add_chrom(i, tpl_ix) for i in range(len(result))])
            window_d.append(window_arr)
    return d, Z, window_d

def get_dstat_df(d, Z, quadruples):
    """
    Takes lists of D values, Z values, and list of name quadruples.
    """
    try:
        dstat_df = pd.DataFrame({'D':d,'Z':Z,
                  'h1':[t[0] for t in quadruples],
                 'h2':[t[1] for t in quadruples],
                 'h3':[t[2] for t in quadruples],
                 'o':[t[3] for t in quadruples]})
    except ValueError, e:
        raise e
        
        
    dstat_df = dstat_df[['h1','h2','h3','o','D','Z']]
    
    return dstat_df

#------------------------------------------------------------------

#------------------------------------------------------------------
#Functions to process D-statistics results
#------------------------------------------------------------------

def dstat_to_pc_mats(dstat_df, maxpar='|D|'):
    pc_df = get_partner_control(dstat_df)
    pc_max_df = get_partner_vs_max_control(pc_df, maxpar)
    dmat, zmat, controlmat = get_partner_vs_max_control_mats(pc_max_df)
    return dmat, zmat, controlmat
        
def get_partner_control(dstat_df):
    dstat_df.loc[:,'|D|'] = dstat_df['D'].apply(abs)
    dstat_df.loc[:,'|Z|'] = dstat_df['Z'].apply(abs)
    #the one that shows geneflow with h3
    dstat_df.loc[:,"p"] = dstat_df["h1"]*(dstat_df["D"]>0) + dstat_df["h2"]*(dstat_df["D"]<0)
    #the other
    dstat_df.loc[:,"c"] = dstat_df["h1"]*(dstat_df["D"]<0) + dstat_df["h2"]*(dstat_df["D"]>0)
    #dstat_df.set_index(["admixture_partner","h3"],inplace=True)
    return dstat_df[['c','p','h3','o','|D|','|Z|']]

def get_partner_vs_max_control(pc_df, maxpar='|D|'):
    """
    maxpar should be |D| or |Z|
    """
    return pc_df.groupby(['p','h3']).apply(lambda df:df.sort(maxpar, ascending=False).iloc[0])

def get_partner_vs_max_control_mats(pc_max_df):
    dmat, zmat, controlmat = pc_max_df['|D|'].unstack(), pc_max_df['|Z|'].unstack(), pc_max_df['c'].unstack()
    return dmat, zmat, controlmat

#------------------------------------------------------------------

#------------------------------------------------------------------
#Functions to plot D-statistics results
#------------------------------------------------------------------



def plot_bubble_chart0(d, z, max_control, order=None, size_factor = 10000, ax=None, transpose=False):
    """
    """
    
    mpl.rcParams['font.size']=16
    
    if ax is None:
        ax = plt.gca()
    stat = '|dstat|'


    if order is not None:
        d = d.loc[order,order]
    
    max_control = max_control.loc[d.index,d.columns]
    z = z.loc[d.index,d.columns]

    if transpose:
        d = d.T
        max_control = max_control.T
        z = z.T

    ax = plt.gca()
    tick_locations = np.arange(0, len(d)+0.5)
    ax.set_xticks(tick_locations)
    ax.set_yticks(tick_locations)
    ax.set_xticklabels(d.columns, rotation='vertical')
    ax.set_yticklabels(d.index)
    ax.set_xlabel(d.columns.name)
    ax.set_ylabel(d.index.name)

    color = []
    x_vals = []
    y_vals = []
    z_scores = []
    sizes = []


    max_z = z.max().max()
    min_z = z.min().min()
    half_z = (max_z - min_z)/2.
    
    pos = 0
    for x in range(max_control.shape[1]):
        for y in range(max_control.shape[0]):
            size = d.iloc[y,x]
            z_score = z.iloc[y,x]
            if not np.isnan(size):
                x_vals.append(x)
                y_vals.append(y)
                sizes.append(size)
                z_scores.append(z_score) 
            if size > 0.001*size_factor/10000:
                try:
                    text = max_control.iloc[y,x][:3]
                    if text == 'ac_':
                        text = max_control.iloc[y,x][3:6]
                    elif text == "A_c":
                        text = max_control.iloc[y,x].split('_')[2][:3]
                    plt.annotate(text,(x,y),horizontalalignment='center',
                                 verticalalignment='center',
                                 color='k' if z_score < half_z else 'w', fontsize=16*np.sqrt(size/0.13*size_factor/10000.))
                except TypeError:
                    pass



    jet = cm = plt.get_cmap('Blues') 
    cNorm  = mpl.colors.Normalize(vmin=min(z_scores), vmax=max(z_scores))

    norm_sizes = np.array(sizes)*size_factor

    bubbles = plt.scatter(np.array(x_vals),np.array(y_vals),
                          c=z_scores,cmap=jet,norm=cNorm,s=norm_sizes)

    cb = plt.colorbar(label="| Z-score |")
    cb.solids.set_rasterized(True) 
    cb.ax.yaxis.labelpad = 10


    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.tick_params(axis=u'both', which=u'both',length=0)
    
    plt.title("Patterson's D {:.2f}% - {:.2f}%".format(d.min().min()*100,d.max().max()*100))
    
    
    return ax



def get_bubble_plot_input(prefix,name,outgroup,tree):
    dstat_df = load_result(prefix, name)
    dstat_df = get_partner_control(dstat_df)
    dstat_df_og = dstat_df[dstat_df['o']==outgroup]
    gene_flow = dstat_df_og.groupby(['admixture_partner','h3']).apply(max_over_control_ingroups)
    pop_tree_order = [c.name for c in tree.get_terminals()]
    order = [o for o in pop_tree_order if o != outgroup]
    return gene_flow, order

def plot_tree_next_to_matrix(tree, outgroup, ax=None):
    if ax is None:
        ax = plt.gca()
    tree_no_outgroup = copy.deepcopy(tree)
    tree_no_outgroup.prune(outgroup)
    draw_tree(tree_no_outgroup,axes=ax0,do_show=False,label_func=lambda s:'')
    
    #ax0.spines['top'].set_visible(False)
    #ax0.spines['right'].set_visible(False)
    #ax0.spines['bottom'].set_visible(False)
    #ax0.spines['left'].set_visible(False)
    #ax0.tick_params(axis=u'both', which=u'both',length=0)
    plt.axis('off')
    
    return ax


#-----------------------------------------------------------------

#------------------------------------------------------------------
#Functions to calculate f-statistics (admixture fraction) and jackknife
#------------------------------------------------------------------

def get_fstat_snpwindow_chrom(chrom, callset, path, populations, 
                              quadruples, controlsamples_h3=1,
                             controlsamples_h2=0, jackknife_window=2000):
    """
    Jackknife window is in number of informative SNPs.
    """
    gen_df = get_gen_df(callset,chrom, path)

    fstats, jackknife_window_fstats = get_fstat_snpwindow(gen_df, quadruples, populations,
                                                         controlsamples_h3=controlsamples_h3,
                                                         controlsamples_h2=controlsamples_h2,
                                                         jackknife_window=2000)
    return fstats, jackknife_window_fstats

def get_fstat_snpwindow(gen_df, quadruples, populations, controlsamples_h3=1,
                        controlsamples_h2=0, jackknife_window=2000):
    
    assert controlsamples_h3>0 or controlsamples_h2>0
    
    fstats = []
    jackknife_window_fstats = []
    
    for j,(h1,h2,h3,o) in enumerate(quadruples):
        #fstat = pd.DataFrame(columns=['num'] \
        #                     + ['denom_h3c_'+str(i) for i in range(controlsamples_h3)] \
        #                     + ['denom_h2c_'+str(i) for i in range(controlsamples_h2)])
        #fstat['num'] = get_numerator(gen_df, (h1,h2,h3,o),
        #                    [populations[p] for p in (h1,h2,h3,o)])
        fstat = get_numerator(gen_df, (h1,h2,h3,o),
                            [populations[p] for p in (h1,h2,h3,o)])
        fstat.name = 'num'
        fstat = pd.DataFrame(fstat)
        for i in range(controlsamples_h3):
            n_samples = len(populations[h3])
            assert n_samples > 1
            s0 = list(np.random.choice(populations[h3], n_samples/2, replace=False))
            s1 = [s for s in populations[h3] if s not in s0]
            denom = get_numerator(gen_df, (h1,h2,h3,o),
                            [populations[h1],s0,s1,populations[o]])
            denom.name = 'denom_h3c_'+str(i)
            fstat = fstat.join(denom, how='outer')
        for i in range(controlsamples_h2):
            #print 'hello'
            n_samples = len(populations[h2])
            assert n_samples > 1
            s0 = list(np.random.choice(populations[h2], n_samples/2, replace=False))
            s1 = [s for s in populations[h2] if s not in s0]
            denom = get_numerator(gen_df, (h1,h2,h3,o),
                            [populations[h1],s0,s1,populations[o]])
            denom.name = 'denom_h2c_'+str(i)
            fstat = fstat.join(denom, how='outer')
        #print fstat
        #attention there seems to be a strange bug in pandas
        #so that df.sum()['a'] != df['a'].sum()
        fstats.append(list(fstat.apply(np.nansum,axis=0)))
        jackknife_window_sum = fstat.rolling(jackknife_window,
                   min_periods=0,
                       center=True).sum().iloc[jackknife_window/2::jackknife_window].dropna()
        jackknife_window_fstats.append(jackknife_window_sum.reset_index(level=1).values.tolist())
    return fstats,jackknife_window_fstats

def get_numerator(gen_df, quadruple, sample_names_quadruple):
    def sample_in_quadruple(sample_name):
        for i in range(4):
            if sample_name in sample_names_quadruple[i]:
                return quadruple[i]
        else:
            return 0
    h1, h2, h3, o = quadruple
    g = gen_df.groupby(sample_in_quadruple,axis=1)
    af = g.mean()/2.
    num = ((af[h1] - af[h2])*(af[h3] - af[o])).dropna()
    num = num[num!=0]
    return num

def reduce_fstat_snpwindow(result, controlsamples_h3, controlsamples_h2):
    
    def fs_from_array(jackknife_arr):
        numdenom_from_jackknife = np.sum(jackknife_arr,axis=0)
        fs_from_jackknife = [numdenom_from_jackknife[0]*1./c for c in numdenom_from_jackknife[1:]]
        return fs_from_jackknife
    
    fs = []
    Zs = []
    #iterate over index of h1,h2,h3,o tuples
    for tpl_ix in range(len(result[0][0])):
        
        numdenom = np.sum(np.array([result[i][0][tpl_ix] for i in range(len(result))]), axis=0)
        fs0 = [numdenom[0]*1./c for c in numdenom[1:]]
        assert len(fs0) == (controlsamples_h3 + controlsamples_h2)
        fs_h3c = fs0[:controlsamples_h3]
        fs_h2c = fs0[controlsamples_h3:]
        
        fs.append([fs_h3c, fs_h2c])


        jackknife_arr = np.concatenate([result[i][1][tpl_ix] for i in range(len(result)) if result[i][1][tpl_ix]])[:,1:]
        
        
        
        fs_from_jackknife = fs_from_array(jackknife_arr)
        

        jackknife_estimates = [fs_from_array(jackknife_arr[np.arange(jackknife_arr.shape[0])!=i]) \
                                                           for i in range(jackknife_arr.shape[0])]
        
        #return fs_from_jackknife, jackknife_estimates
        print "Number of jackknife estimates for quadruple {}:".format(tpl_ix), len(jackknife_estimates)
        
        zscores = fs_from_jackknife/(np.std(jackknife_estimates,axis=0,ddof=1)*np.sqrt(len(jackknife_estimates)))
        Zs.append([list(zscores[:controlsamples_h3]),list(zscores[controlsamples_h3:])])

    return fs, Zs


