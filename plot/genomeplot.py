import matplotlib as mpl
from matplotlib import pyplot as plt
import numpy as np
import pandas

def plot_along_chrom(df,chrom,quantile,window_size):
    chr_len = chrom_length.ix[chrom]
    n_subplot = len(pops)
    #print round(window_size/1000)
    fig = plt.figure(1, figsize=(15./(1-log(10000./window_size))\
                           *chr_len/chrom_length.ix["CAE1"],1.5*n_subplot))
    plt.suptitle("Akey's d, "+str(int(window_size/1000))+"kb windows",fontsize="15") 
    hspace = 0.3
    fig.subplots_adjust(hspace=hspace)
    for i,pop in enumerate(pops):
        qt = quantile[pop]
        ax = fig.add_subplot(n_subplot,1,i+1)
        ax.set_frame_on(False)
        ax.set_yticks([0,round(qt)])
        if i < len(pops)-1:
            ax.set_xticks([])
            ax.get_xaxis().set_visible(False)
        else:
            ax.set_xlabel(chrom)
        #total_3window_d[idx][pop].plot(style=".",markersize=5,color=colors[i],ax=ax)
        ax.plot(window_d.ix[chrom][pop].index,window_d.ix[chrom][pop].values,".",markersize=5,color=colors[i])
        ax.set_title(pop,position=(0.5,0.95))
        ax.set_xlim([0,chr_len+1000])
    return fig

def get_chrom_grid(chrom_lens=None,parent_gridelement=None):
    """
    subplot_spec is the parent gird
    """
    if chrom_lens is None:
        #only autosomes
        chrom_lens = np.array([126035930,  90373283,  92142175,  91010382,  75399963,  50890351,
       135778131, 139301422, 125710982, 128595539, 128539186, 108555830,
        98384682, 107702431,  91754291,  75148670,  71996105,  72318688,
        33263144, 130588469, 127223203, 101219884,  82825804,  84932903,
        85787240,  58131712,  48547382,  21531802,  24206276, 130038232,
         6181219])[:-2]
    if parent_gridelement is None:
        parent_gridelement = mpl.gridspec.GridSpec(1, 1)[0] 
    gs = mpl.gridspec.GridSpecFromSubplotSpec(1, len(chrom_lens),
                                                     width_ratios=chrom_lens.astype(float)/sum(chrom_lens),
                                                                                 subplot_spec=parent_gridelement)
    return gs, parent_gridelement

def get_chrom_axes(grid=None,parent_gridelement=None,fig=None,
                      ylabel = "", plot_xlabels=True, color=None):
    """
    axes of the returned figure are the chromosomes 
    plus one global axes across all chromosomes
    """
    if plotfun is None:
        plotfun = plt.plot
    if grid is None:
        grid, parent_gridelement = get_chrom_grid(parent_gridelement=parent_gridelement)
    if fig is None:
        fig = plt.figure(1,figsize=(20,2))
    fig.subplots_adjust(wspace=0)
    #print color
    if color is None:
        print "changed to blue"
        color = mpl.rcParams['axes.color_cycle'][0]
    color = converter.to_rgb(color) 
    chroms = series.index.get_level_values(level=0).unique()
    for i,(chrom, gr) in enumerate(zip(chroms,grid)):
        if i == 0:
            ax = fig.add_subplot(gr)
        else:
            ax = fig.add_subplot(gr, sharey=ax)
            ax.set_yticks([])
        
        if i % 2 == 0:
            color1 = tuple(np.array(color)*0.6)
            if plot_xlabels:
                ax.set_xlabel(chrom.split("E")[1])
        else:
            color1 = color
        #print color1
        ax.set_frame_on(False)
        ax.set_xticks([])
        ax.set_color_cycle([color1])
    ylim = ax.get_ylim()
    #we could also get the parent element with gr.get_topmost_subplotspec()
    total_ax = fig.add_subplot(parent_gridelement)
    total_ax.set_frame_on(False)
    total_ax.spines['top'].set_visible(True)
    total_ax.set_ylim(ylim)
    total_ax.set_xticks([])
    if not ylabel and series.name is not None:
        ylabel = series.name
    total_ax.set_ylabel(ylabel,fontsize = 14)
    if plot_xlabels:
        total_ax.set_xlabel("Chromosomes",labelpad=25,fontsize = 14)
    return fig.get_axes()

converter =  mpl.colors.ColorConverter()
def plot_chrom_series(series,chrom_len,plotfun=None,grid=None,parent_gridelement=None,fig=None,
                      ylabel = "", plot_xlabels=True, color=None,title=None,rightlabel=None,**kwa):
    """

    chrom_len ... series or dict with chrom names as keys and chromosomes length as values
    axes of the returned figure are the chromosomes 
    plus one global axes across all chromosomes

    
    """
    if plotfun is None:
        plotfun = plt.plot
    if grid is None:
        grid, parent_gridelement = get_chrom_grid(chrom_lens=chrom_len,parent_gridelement=parent_gridelement)
    if fig is None:
        fig = plt.figure(1,figsize=(20,2))
    fig.subplots_adjust(wspace=0)
    #print color
    if color is None:
        color = mpl.rcParams['axes.color_cycle'][0]
    color = converter.to_rgb(color) 
    #chroms = series.index.get_level_values(level=0).unique()
    try: 
        chroms = chrom_len.index
    except AttributeError:
        chroms = chrom_len.keys()

    axes = {}
    for i,(chrom, gr) in enumerate(zip(chroms,grid)):
        #print chrom
        if i == 0:
            ax = fig.add_subplot(gr)
        else:
            ax = fig.add_subplot(gr, sharey=ax)
            ax.set_yticks([])

        axes.update({chrom:ax})
        
        if i % 2 == 0:
            color1 = tuple(np.array(color)*0.6)
            if plot_xlabels:
                ax.set_xlabel(chrom)
        else:
            color1 = color
        #print color1
        ax.set_frame_on(False)
        ax.set_xticks([])
        ax.set_xlim([0,chrom_len[chrom]])
       # if chrom == 'CAE3':
            #print series.ix[chrom].index.values[~np.isnan(series.ix[chrom].values)]
        try:
            chrom_series = series.loc[chrom]
        except KeyError:
            continue
        plotfun(chrom_series.index.values,chrom_series.values,'.',color=color1,rasterized=True,**kwa)
    ylim = ax.get_ylim()
    #we could also get the parent element with gr.get_topmost_subplotspec()
    #we don't share y here because otherwise the ticks are interdependent,
    #is there a way to turn tick on for only one of the shared axes?
    total_ax = fig.add_subplot(parent_gridelement) 
    axes.update({"total_ax" : total_ax})
    total_ax.set_frame_on(False)
    total_ax.spines['top'].set_visible(True)
    total_ax.set_ylim(ylim)
    total_ax.set_xticks([])
    if not ylabel and series.name is not None:
        ylabel = series.name
    total_ax.set_ylabel(ylabel,fontsize = 14)
    total_ax.set_yticks
    
    if plot_xlabels:
        total_ax.set_xlabel("Chromosomes",labelpad=25,fontsize = 14)
        
    if rightlabel is not None:
        ax2 = ax.twinx()
        axes.update({"right_ax":ax2})
        ax2.set_ylabel(rightlabel, color=color1,fontsize=16)
        ax2.set_yticks([])
        ax2.spines['top'].set_visible(False)
        ax2.spines['right'].set_visible(False)
        ax2.spines['bottom'].set_visible(False)
        ax2.spines['left'].set_visible(False)
        
    if title is not None:
        plt.title(title)#, position=(0.5,1.02),va='bottom'

    return fig, axes, grid, parent_gridelement

def plot_features(gene_df,annotate=True,feature_name="symbol",ax=None,xlim=None):
    if ax is None:
        ax = plt.gca()
        
    gene_df.sort('start', inplace=True)
    bar_thickness = 0.05
    text_space = 0.05
    text_dist = 0.01
    row_height = bar_thickness + text_space
    
    row = 0
    tot_rows = 0
    max_rows = 3
    prev_end = -np.inf
    for i, g in gene_df.iterrows():
        if g['start'] <= prev_end and row < max_rows:
            row += 1
        else:
            row = 0
        if annotate:
            ax.annotate(g[feature_name], ((g['start']+g['end'])/2.,-bar_thickness-row*row_height-text_dist), #
                          size=16, ha="right", va="top",rotation=30,color='k')#rotation = 270,+g['length']/2.
        plt.barh(bottom= -bar_thickness-row*row_height,
                 left=g['start'],
                 width=g['end']-g['start'],
                 height=bar_thickness,
                 color = 'yellow', alpha=1) #0.5
        tot_rows = max(tot_rows, row)
        prev_end = g['end']
        
    
           
  
    ax.set_ylim([-(tot_rows+1)*row_height,0])
    
    if xlim is None:
        ax.set_xlim([gene_df['start'].min(),gene_df['end'].max()])
    else:
        ax.set_xlim(xlim)
    ax.set_axis_off()
    
    return ax
