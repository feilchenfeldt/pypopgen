from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio.Phylo.TreeConstruction import _DistanceMatrix

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
