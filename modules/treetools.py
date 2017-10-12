import StringIO
import ete3
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio.Phylo.TreeConstruction import _DistanceMatrix

import matplotlib as mpl

#class Tree(self):

#    def __init__(self,):


def dm_to_tree(dm, ladderize=True):
    dm = dm.astype(float)
    distance_triangular = [list(dm.values[i,:i+1]) for i in range(len(dm))]
    try:
        dm = _DistanceMatrix(names= [str(i) for i in dm.columns],
                    matrix=distance_triangular)
    except Exception,e:
        print list(dm.columns)
        print [type(i) for i in dm.columns]
        print type(distance_triangular)
        print  type(distance_triangular[0])
        print  set([str(type(i)) for j in distance_triangular for i in j])
        print distance_triangular
        raise e
    constructor = DistanceTreeConstructor()
    tree = constructor.nj(dm)
    for c in tree.get_nonterminals():
        c.name = None
    if ladderize:
        tree.ladderize()
    return tree

def get_local_tree(chrom, start, end, vcf_fn, samples=None, outgroup=None, plot=False):
    pwd = hap.get_pairwise_diff(vcf_fn, chrom=chrom, start=start, end=end, samples=samples, chunksize=30000)
    distance_triangular = [list(pwd.values[i,:i+1]) for i in range(len(pwd))]
    dm = _DistanceMatrix(names= list(pwd.columns),
                        matrix=distance_triangular)
    constructor = DistanceTreeConstructor()
    tree = constructor.nj(dm)
    if outgroup is not None:
        tree.root_with_outgroup(outgroup)
    tree.ladderize()
    for t in tree.get_nonterminals():
        t.name = None
    tree_no = copy.deepcopy(tree)
    if outgroup is not None:
        tree_no.prune(outgroup)
    #tree_no.prune('AstTwe1')
    #tree_no.prune('SerRob1')
    if plot:
        fig = plt.figure(figsize=(15,50))
        ax = plt.gca()

        Phylo.draw(tree_no, axes=ax, do_show=False)
        ax.spines['left'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.set_yticks([])
        ax.set_ylabel('')
    return pwd, tree

def phylo_to_ete(phylo_tree):
    treeio = StringIO.StringIO()
    Phylo.write(phylo_tree, treeio, format='newick')
    ete_tree = ete3.Tree(treeio.getvalue())
    #make sure no additional outer quotes in names with whitespace
    for t in ete_tree.get_leaves():
            t.name = t.name.strip("' ")
    return ete_tree


def phylo_from_str(tree_str):
    treeio = StringIO.StringIO(tree_str)
    phylo_tree = Phylo.read(treeio, format='newick')
    return phylo_tree

def consistent_with_tree(ete_tree, h1, h2, h3, h4):
    """"
    Returns True if the lineages have relationship
       /-h4  
      /
     |   /-h3
     |--|
         |   /-h2
          \-|
             \-h1

    and h1 != h2.
    Returns False otherwise.
    """
    if h1 == h2:
        return False
    try:
        h1h2 = ete_tree.get_common_ancestor([h1,h2])
        h1h2h3 = ete_tree.get_common_ancestor([h1,h2,h3])
        h1h2h3h4 = ete_tree.get_common_ancestor([h1,h2,h3,h4])
    except ValueError:
        #One of the samples not in tree.
        return False
    return bool(h1h2.get_distance(h1h2h3, topology_only=True)) \
            & bool(h1h2h3.get_distance(h1h2h3h4, topology_only=True))

def get_consistent_df(stat_df, ete_tree):
        """
        Get a data frame with the subset
        of tuples that are consistent with a
        given ete tree.
        Parameters:
        stat_df : data frame with f4-like statistics
                  where index levels are [h1, h2, h3, h4]
        ete_tree : ete3 tree object of all samples.
                   Needs to be rooted and include
                   all outgroups.
        
        """
        consistent_tpls = [consistent_with_tree(ete_tree, *t) for t in stat_df.index.values]
        return stat_df[consistent_tpls]

#-------------------------------------------
###TREE PLOTTING AND VISUALISATION#########
#-------------------------------------------

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
                               heights[clade.clades[-1]]) / 2.0)  # -1 #hs subtract -1

        if tree.root.clades:
            calc_row(tree.root)
        return {k: v - 1 for (k, v) in heights.iteritems()}  # HS -1

    x_posns = get_x_positions(tree)
    y_posns = get_y_positions(tree)

    # The function draw_clade closes over the axes object
    if axes is None:
        fig = mpl.pyplot.figure()
        axes = fig.add_subplot(1, 1, 1)
    elif not isinstance(axes, mpl.pyplot.matplotlib.axes.Axes):
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
            lw = clade.width * mpl.pyplot.rcParams['lines.linewidth']
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
            y_top = y_posns[clade.clades[0]]  # - 1
            y_bot = y_posns[clade.clades[-1]]  # - 1
            # Only apply widths to horizontal lines, like Archaeopteryx
            draw_clade_lines(use_linecollection=True, orientation='vertical',
                             x_here=x_here, y_bot=y_bot, y_top=y_top, color=color, lw=lw)
            # Draw descendents
            for child in clade:
                draw_clade(child, x_here, color, lw)

    draw_clade(tree.root, 0, 'k', mpl.pyplot.rcParams['lines.linewidth'])

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
    axes.set_ylim(0, max(y_posns.values()))

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
            getattr(mpl.pyplot. str(key))(**dict(value))
        elif not (isinstance(value[0], tuple)):
            getattr(mpl.pyplot. str(key))(*value)
        elif (isinstance(value[0], tuple)):
            getattr(mpl.pyplot. str(key))(*value[0], **dict(value[1]))
    return axes

def plot_node_tree(tree, ax=None, x0=0,y0=0, plot_leaf_names=True, em=0.7,fontsize=12):
    if ax is None:
        ax = mpl.pyplot.gca()
    x = x0
    y = y0
    maxletters = max([len(n) for n in tree.get_leaf_names()])
    if not plot_leaf_names:
        em = 0


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
                ax.annotate(node.name,xy=(depth+em*maxletters,y),fontsize=12, va='center',ha='right')
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

    mpl.pyplot.axis('off')

    return ax



