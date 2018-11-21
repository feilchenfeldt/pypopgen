"""
Some functions to plot phylogenetic
trees.
hannes.svardal@gmail.com
"""
import copy
import matplotlib as mpl
import numpy as np
from matplotlib import pyplot as plt


def plot_node_tree(tree, ax=None, x0=0,y0=0, plot_leaf_names=True, em=0.7,fontsize=12):
    """
    Plot a special tree that has
    a dashed line going from each
    internal node to the present.
    This is used as a legend to visualise
    data with values for both internal
    and external nodes. 
    """
    if ax is None:
        ax = plt.gca()
    x = x0
    y = y0
    maxletters = max([len(n) for n in tree.get_leaf_names()])
    if not plot_leaf_names:
        em = 0
    
    visited = []
    depth = tree.get_farthest_leaf(topology_only=True)[1] +2
    for node in tree.iter_descendants('preorder'):
        visited.append(node)

        x = node.get_distance(node.get_tree_root(),topology_only=True)


        hline, = ax.plot([x,x+1], [y,y], '-',color='k')

        try:
            rows_below = len(node.get_descendants())
        except IndexError:
            rows_below = 0

        if node.is_leaf():


            nletters = len(node.name)
            if plot_leaf_names:
                ax.annotate(node.name,xy=(depth+em*maxletters,y), va='center',ha='right', fontsize=fontsize)
            hline2, = ax.plot([x+1,depth+em*(maxletters-nletters)], [y,y], '-',color='k')
        else:
            #print rows_below
            try:
                sister = node.get_sisters()[0]
                if sister not in visited:
                    vline, = ax.plot([x,x], [y,y-rows_below-1], '-',color='k')
            except IndexError:
                pass

            hline2, = ax.plot([x+1,depth++em*maxletters], [y,y], ls='dotted',color='k')
            vline2, = ax.plot([x+1,x+1], [y,y-2], '-',color='k')

        y = y - 1

    plt.axis('off')
    
    return ax

def draw_tree(tree, x , y, depth,ax=None, orientation='horizontal',
             x_increment=1, y_increment=1,
             topology_only=True, plot_leaf_names=True,
             name_margin=7, correction=0.2):
    """
    draw a phylogenetic tree in ete3 format
    
    currently only plots topology,
    implementing distance should be easy
    
    topology_only=False not implemented yet
    
    
    
    """    
    def get_x_y(x,y,orientation):
        if orientation == 'horizontal':
            return x,y
        elif orientation == 'vertical':
            return y,x
        
    x_inc, y_inc = get_x_y(x_increment, y_increment, orientation)
    
    
    if ax is None:
        ax = plt.gca()
    
    if not tree.is_leaf():
        c0, c1 = tree.get_children()
        #node_height = np.ceil(len(c0.get_leaves())/2.) + np.ceil(len(c1.get_leaves())/2.)
        try:
            node_height0 = len(c0.get_children()[1].get_leaves())
        except IndexError:
            node_height0 = 0.5
        try:
            node_height1 = len(c1.get_children()[0].get_leaves())
        except IndexError:
            node_height1 = 0.5
        
        hline, = ax.plot(*get_x_y([x,x+x_inc],[y,y],orientation),color='k')

        vline, = ax.plot(*get_x_y([x+x_inc,x+x_inc], [y,y+(node_height0)*y_inc],orientation),
                                                                         color='k')
        draw_tree(c0, x+x_inc, y+(node_height0)*y_inc, depth=depth, ax=ax,
                    orientation=orientation, x_increment=x_increment,y_increment=y_increment,
                  name_margin=name_margin,correction=correction)

        vline, = ax.plot(*get_x_y([x+x_inc,x+x_inc], [y,y-(node_height1)*y_inc],orientation),
                                                                     color='k')
        draw_tree(c1, x+x_inc, y-(node_height1)*y_inc, depth=depth, ax=ax, 
                    orientation=orientation, x_increment=x_increment,y_increment=y_increment,
                             name_margin=name_margin,correction=correction)

    else:
        
        if plot_leaf_names:
            if orientation=='vertical':
                if x_inc>0:
                    rotation = 270
                    va = 'top'
                else:
                    rotation = 90
                    va = 'bottom'
                ha = 'center'
            else:
                rotation = 0
                va ='center'
                if x_inc>0:
                    ha = 'right'
                else:
                    ha = 'left'
            
            
            lbl = ax.annotate(tree.name,xy=get_x_y((depth+name_margin)*x_inc,y,orientation),
                        fontsize=10, va=va,ha=ha,rotation=rotation)
            r = ax.get_figure().canvas.get_renderer()
            tbox = lbl.get_window_extent(renderer=r)
            dbox = tbox.transformed(ax.transData.inverted())
            text_width = dbox.x1-dbox.x0
            text_height = dbox.y1-dbox.y0
            text_span = get_x_y(text_width ,text_height ,orientation)[0]

            
            hline, = ax.plot(*get_x_y([x,(depth+name_margin-text_span+correction)*x_inc],
                                      [y,y],orientation),
                                                                                 color='k')  
    return ax


def plot_tree(tree0, x0=0, y0=0, plot_labels=True, ax=None, labeldist_correct_factor=1,fontsize=10): 
    tree = copy.deepcopy(tree0)
    depth = tree.get_farthest_leaf(topology_only=True)[1] +2
    use_edge_length = False
    xstep = 1
    ystep = 1
    
    terminals = tree.get_leaves()
    
    if ax is None:
        ax = plt.gca()
    
    
    ylim_min = y0

    ylim_max = y0 + depth + 9

    label_ymaxs = []

    xcoords = np.arange(x0,x0+len(terminals)*xstep,xstep)
    if plot_labels:
        anns = []
        for x, t in zip(xcoords, 
                       terminals):
            anns.append(ax.annotate(t.name, xy=(x, y0), rotation=90, 
                        horizontalalignment='center', verticalalignment='bottom',
                        fontsize=fontsize
                                   # xycoords='figure pixels'
                                   ))

        #ax.annotate

        ax.figure.canvas.draw()

        for ann in anns:
            bbox =  ann.get_window_extent()
            (xmin,ymin),(xmax,ymax) = ax.transData.inverted().transform(bbox)
            #print bbox
            #print xmin, xmax, ymin, ymax
            
            label_ymaxs.append(y0+(ymax-ymin)*(ylim_max-ylim_min))
            
        ylim_max_old = ylim_max
        ylim_max = depth*ystep + max(label_ymaxs) + ystep
        for i in range(len(label_ymaxs)):
            label_ymaxs[i] = (label_ymaxs[i] + label_ymaxs[i]*(ylim_max - ylim_max_old)/ylim_max_old)*labeldist_correct_factor
            
    else:
        label_ymaxs = [y0]*len(terminals)


    maxlabelpos = max(label_ymaxs)

    
    
    for x, yl, t in zip(xcoords, label_ymaxs, terminals):
        t.add_feature('x',x)
        t.add_feature('y',maxlabelpos)
        ax.plot([x,x],[yl, maxlabelpos], color='k')
        #height = t.up.get_farthest_leaf(topology_only=True)[1]+1
    
    while tree:
        for t in tree.get_leaves():
            sisters = t.get_sisters()
            assert len(sisters) <= 1, t
            if not sisters:
                #print "no sisters:", t
                continue
            s = sisters[0]
            if s.is_leaf():
                #print "plotting", t.name, s.name
                ypos = max(t.y, s.y) + ystep
                ax.plot([t.x,t.x],[t.y, ypos], color='k')
                ax.plot([s.x,s.x],[s.y, ypos], color='k')
                ax.plot([t.x,s.x],[ypos, ypos], color='k')
                #p = t.up
                #p.add_feature('x',(t.x+s.x)/2.)
                #p.add_feature('y',ypos)
                t.x = (t.x+s.x)/2.
                t.y = ypos
                s.delete()
                
        #print len(tree), tree
        if len(tree)<2:
            break
            
    



    ax.set_ylim([ylim_min,ylim_max])
    ax.set_xlim([x0-0.5, x0+len(terminals)*xstep+0.5])
    
    #plt.axis('off')
    #fig.axes.get_xaxis().set_visible(False)
    #fig.axes.get_yaxis().set_visible(False)
    
    return ax
