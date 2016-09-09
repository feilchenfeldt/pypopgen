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
                   min_periods=0,#what should this be?
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
    draw_tree(tree_no_outgroup,axes=ax,do_show=False,label_func=lambda s:'')
    
    #ax0.spines['top'].set_visible(False)
    #ax0.spines['right'].set_visible(False)
    #ax0.spines['bottom'].set_visible(False)
    #ax0.spines['left'].set_visible(False)
    #ax0.tick_params(axis=u'both', which=u'both',length=0)
    plt.axis('off')
    
    return ax

def draw_tree(tree, label_func=str, do_show=True, show_confidence=True,
         # For power users
         axes=None, branch_labels=None, *args, **kwargs):
    """ HS: This function is taken from Bio.Phylo._utils
        and modified to return axis and to plot 
    
    Plot the given tree using matplotlib (or pylab).

    The graphic is a rooted tree, drawn with roughly the same algorithm as
    draw_ascii.

    Additional keyword arguments passed into this function are used as pyplot
    options. The input format should be in the form of:
    pyplot_option_name=(tuple), pyplot_option_name=(tuple, dict), or
    pyplot_option_name=(dict).

    Example using the pyplot options 'axhspan' and 'axvline':

    >>> Phylo.draw(tree, axhspan=((0.25, 7.75), {'facecolor':'0.5'}),
    ...     axvline={'x':'0', 'ymin':'0', 'ymax':'1'})

    Visual aspects of the plot can also be modified using pyplot's own functions
    and objects (via pylab or matplotlib). In particular, the pyplot.rcParams
    object can be used to scale the font size (rcParams["font.size"]) and line
    width (rcParams["lines.linewidth"]).

    :Parameters:
        label_func : callable
            A function to extract a label from a node. By default this is str(),
            but you can use a different function to select another string
            associated with each node. If this function returns None for a node,
            no label will be shown for that node.
        do_show : bool
            Whether to show() the plot automatically.
        show_confidence : bool
            Whether to display confidence values, if present on the tree.
        axes : matplotlib/pylab axes
            If a valid matplotlib.axes.Axes instance, the phylogram is plotted
            in that Axes. By default (None), a new figure is created.
        branch_labels : dict or callable
            A mapping of each clade to the label that will be shown along the
            branch leading to it. By default this is the confidence value(s) of
            the clade, taken from the ``confidence`` attribute, and can be
            easily toggled off with this function's ``show_confidence`` option.
            But if you would like to alter the formatting of confidence values,
            or label the branches with something other than confidence, then use
            this option.
    """
    import matplotlib.collections as mpcollections

    # Arrays that store lines for the plot of clades
    horizontal_linecollections = []
    vertical_linecollections = []

    # Options for displaying branch labels / confidence
    def conf2str(conf):
        if int(conf) == conf:
            return str(int(conf))
        return str(conf)
    if not branch_labels:
        if show_confidence:
            def format_branch_label(clade):
                if hasattr(clade, 'confidences'):
                    # phyloXML supports multiple confidences
                    return '/'.join(conf2str(cnf.value)
                                    for cnf in clade.confidences)
                if clade.confidence:
                    return conf2str(clade.confidence)
                return None
        else:
            def format_branch_label(clade):
                return None
    elif isinstance(branch_labels, dict):
        def format_branch_label(clade):
            return branch_labels.get(clade)
    else:
        assert callable(branch_labels), \
            "branch_labels must be either a dict or a callable (function)"
        format_branch_label = branch_labels

    # Layout

    def get_x_positions(tree):
        """Create a mapping of each clade to its horizontal position.

        Dict of {clade: x-coord}
        """
        depths = tree.depths()
        # If there are no branch lengths, assume unit branch lengths
        if not max(depths.values()):
            depths = tree.depths(unit_branch_lengths=True)
        return depths

    def get_y_positions(tree):
        """Create a mapping of each clade to its vertical position.

        Dict of {clade: y-coord}.
        Coordinates are negative, and integers for tips.
        """
        maxheight = tree.count_terminals()
        # Rows are defined by the tips
        heights = dict((tip, maxheight - i)
                       for i, tip in enumerate(reversed(tree.get_terminals())))

        # Internal nodes: place at midpoint of children
        def calc_row(clade):
            for subclade in clade:
                if subclade not in heights:
                    calc_row(subclade)
            # Closure over heights
            heights[clade] = ((heights[clade.clades[0]] +
                              heights[clade.clades[-1]]) / 2.0) #-1 #hs subtract -1
            
        if tree.root.clades:
            calc_row(tree.root)
        return {k:v-1 for (k,v) in heights.iteritems()} #HS -1

    x_posns = get_x_positions(tree)
    y_posns = get_y_positions(tree)
    
    # The function draw_clade closes over the axes object
    if axes is None:
        fig = plt.figure()
        axes = fig.add_subplot(1, 1, 1)
    elif not isinstance(axes, plt.matplotlib.axes.Axes):
        raise ValueError("Invalid argument for axes: %s" % axes)

    def draw_clade_lines(use_linecollection=False, orientation='horizontal',
                         y_here=0, x_start=0, x_here=0, y_bot=0, y_top=0,
                         color='black', lw='.1'):
        """Create a line with or without a line collection object.

        Graphical formatting of the lines representing clades in the plot can be
        customized by altering this function.
        """
        if (use_linecollection is False and orientation == 'horizontal'):
            axes.hlines(y_here, x_start, x_here, color=color, lw=lw)
        elif (use_linecollection is True and orientation == 'horizontal'):
            horizontal_linecollections.append(mpcollections.LineCollection(
                [[(x_start, y_here), (x_here, y_here)]], color=color, lw=lw),)
        elif (use_linecollection is False and orientation == 'vertical'):
            axes.vlines(x_here, y_bot, y_top, color=color)
        elif (use_linecollection is True and orientation == 'vertical'):
            vertical_linecollections.append(mpcollections.LineCollection(
                [[(x_here, y_bot), (x_here, y_top)]], color=color, lw=lw),)

    def draw_clade(clade, x_start, color, lw):
        """Recursively draw a tree, down from the given clade."""
        x_here = x_posns[clade]
        y_here = y_posns[clade]
        # phyloXML-only graphics annotations
        if hasattr(clade, 'color') and clade.color is not None:
            color = clade.color.to_hex()
        if hasattr(clade, 'width') and clade.width is not None:
            lw = clade.width * plt.rcParams['lines.linewidth']
        # Draw a horizontal line from start to here
        draw_clade_lines(use_linecollection=True, orientation='horizontal',
                         y_here=y_here, x_start=x_start, x_here=x_here, color=color, lw=lw)
        # Add node/taxon labels
        label = label_func(clade)
        if label not in (None, clade.__class__.__name__):
            axes.text(x_here, y_here, ' %s' %
                      label, verticalalignment='center')
        # Add label above the branch (optional)
        conf_label = format_branch_label(clade)
        if conf_label:
            axes.text(0.5 * (x_start + x_here), y_here, conf_label,
                      fontsize='small', horizontalalignment='center')
        if clade.clades:
            # Draw a vertical line connecting all children
            y_top = y_posns[clade.clades[0]] #- 1
            y_bot = y_posns[clade.clades[-1]] #- 1
            # Only apply widths to horizontal lines, like Archaeopteryx
            draw_clade_lines(use_linecollection=True, orientation='vertical',
                             x_here=x_here, y_bot=y_bot, y_top=y_top, color=color, lw=lw)
            # Draw descendents
            for child in clade:
                draw_clade(child, x_here, color, lw)
                
                
    draw_clade(tree.root, 0, 'k', plt.rcParams['lines.linewidth'])

    # If line collections were used to create clade lines, here they are added
    # to the pyplot plot.
    for i in horizontal_linecollections:
        axes.add_collection(i)
    for i in vertical_linecollections:
        axes.add_collection(i)

    # Aesthetics

    if hasattr(tree, 'name') and tree.name:
        axes.set_title(tree.name)
    axes.set_xlabel('branch length')
    axes.set_ylabel('taxa')
    # Add margins around the tree to prevent overlapping the axes
    xmax = max(x_posns.values())
    axes.set_xlim(-0.05 * xmax, 1.25 * xmax)
    # Also invert the y-axis (origin at the top)
    # Add a small vertical margin, but avoid including 0 and N+1 on the y axis
    #axes.set_ylim(max(y_posns.values()) + 0.8, 0.2)
    axes.set_ylim(0,max(y_posns.values()))

    # Parse and process key word arguments as pyplot options
    for key, value in kwargs.items():
        try:
            # Check that the pyplot option input is iterable, as required
            [i for i in value]
        except TypeError:
            raise ValueError('Keyword argument "%s=%s" is not in the format '
                             'pyplot_option_name=(tuple), pyplot_option_name=(tuple, dict),'
                             ' or pyplot_option_name=(dict) '
                             % (key, value))
        if isinstance(value, dict):
            getattr(plt, str(key))(**dict(value))
        elif not (isinstance(value[0], tuple)):
            getattr(plt, str(key))(*value)
        elif (isinstance(value[0], tuple)):
            getattr(plt, str(key))(*value[0], **dict(value[1]))
    return axes   



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
                                                         jackknife_window=jackknife_window)
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
            assert n_samples > 1, "cannot split p3 {}:{}".format(h3,populations[h3])
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

def get_fstat_df(fs, Zs, quadruples):
    """
    Takes lists of lists of f values, Z values, and list of name quadruples.
    """
    d = {'h1':[t[0] for t in quadruples],
                 'h2':[t[1] for t in quadruples],
                 'h3':[t[2] for t in quadruples],
                 'o':[t[3] for t in quadruples]}
    d.update({'f_h3c_'+str(i):[fs[t][0][i] for t in range(len(quadruples))] for i in range(len(fs[0][0]))})
    
    d.update({'Z_h3c_'+str(i):[Zs[t][0][i] for t in range(len(quadruples))] for i in range(len(fs[0][0]))})
    
    d.update({'f_h2c_'+str(i):[fs[t][1][i] for t in range(len(quadruples))] for i in range(len(fs[0][1]))})
    d.update({'Z_h2c_'+str(i):[Zs[t][1][i] for t in range(len(quadruples))] for i in range(len(fs[0][1]))})
    try:
        fstat_df = pd.DataFrame(d)
    except ValueError, e:
        raise e
    
    return fstat_df

def get_partner_control_f(fstat_df):
    """
    assumes a column f which has to be computed before
    also assumes a column Z
    """
    fstat_df = fstat_df.copy()
    fstat_df.loc[:,'|f|'] = fstat_df['f'].apply(abs)
    fstat_df.loc[:,'|Z|'] = fstat_df['Z'].apply(abs)
    #the one that shows geneflow with h3
    fstat_df.loc[:,"c"] = fstat_df["h1"]*(fstat_df["f"]>0) + fstat_df["h2"]*(fstat_df["f"]<0)
    #the other
    fstat_df.loc[:,"p"] = fstat_df["h1"]*(fstat_df["f"]<0) + fstat_df["h2"]*(fstat_df["f"]>0)
    return fstat_df[['c','p','h3','o','|f|','|Z|']]

#------------------------------------------------------------------
#Functions to plot f-statistics results
#------------------------------------------------------------------





#------------------------------------------------------------------
#F-reduced
#------------------------------------------------------------------

def get_f_reduced(f_df_pc, etetree, outgroup):
    t = copy.deepcopy(etetree)
    f_reduced = f_df_pc.copy()
    for node in  t.iter_descendants():
        if node.children:
            l = node.children[0]
            r = node.children[1]
            lleaves = l.get_leaf_names()
            rleaves = r.get_leaf_names()
            for h3 in [p for p in pops if p not in lleaves+rleaves and p!=outgroup]:
                for p in lleaves:
                    for c in rleaves:
                        control_df = f_df_pc[(f_df_pc['p']==p)&(f_df_pc['c'].isin(lleaves))&(f_df_pc['h3']==h3)]
                        if len(control_df):
                            control = control_df.sort_values('|f|',ascending=False).iloc[0]
                            f_reduced.loc[(f_df_pc['p']==p)&(f_df_pc['c']==c)&(f_df_pc['h3']==h3),'|f|'] = \
                            f_reduced.loc[(f_df_pc['p']==p)&(f_df_pc['c']==c)&(f_df_pc['h3']==h3),'|f|'] -control['|f|']
                for p in rleaves:
                    for c in lleaves:
                        control_df = f_df_pc[(f_df_pc['p']==p)&(f_df_pc['c'].isin(rleaves))&(f_df_pc['h3']==h3)]
                        if len(control_df):
                            control = control_df.sort_values('|f|',ascending=False).iloc[0]
                            f_reduced.loc[(f_df_pc['p']==p)&(f_df_pc['c']==c)&(f_df_pc['h3']==h3),'|f|'] = \
                            f_reduced.loc[(f_df_pc['p']==p)&(f_df_pc['c']==c)&(f_df_pc['h3']==h3),'|f|'] -control['|f|']
    return f_reduced, t


