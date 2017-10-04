import StringIO
import ete3
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio.Phylo.TreeConstruction import _DistanceMatrix

#class Tree(self):

#    def __init__(self,):


def dm_to_tree(dm):
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
