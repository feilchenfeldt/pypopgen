import numpy as np
import pandas as pd
import scipy
import treetools as tt
import warnings


def get_ne(pwd, mutation_rate, accessible_size, groupings=None):
    """
    If groupings are given. The Ne represents 
    the Ne of the (stuctured) populations 
    formed by the groupings.

    Returns:
    Dictionary of group names and Ne.
    """
    assert (pwd.index == pwd.columns).all()
    if groupings is None:
        pwd_grouped = pwd
    else:
        pwd_grouped = pwd.groupby(groupings).mean().groupby(groupings,axis=1).mean()
    
    return {k:v for k,v in zip(pwd_grouped.index,
                               np.diagonal(pwd_grouped/(4.*mutation_rate*accessible_size)))}

def get_group_pwd(pwd, groupings):
    """
    Note the difference to get_ne.
    Get ne averages across comparisons
    of individuals within groups.
    This averages over cross-group
    comparisons.

    Returns: DataFrame of pairwise differences
    between groups.
    """
    return pwd.groupby(groupings)\
                        .apply(lambda df:df.groupby(groupings,axis=1)\
                               .apply(np.diagonal)).applymap(np.mean)

def get_samples_per_group(pwd, groupings=None, haploid=True):
    """
    Returns the number of samples in each group.
    """
    if groupings is None:
        return {n:(1+haploid) for n in pwd.index}
    else:
        return (pwd.groupby(groupings).apply(len)*(1+haploid)).to_dict()

def get_split_diff(pwd):
    """
    Returns an estimate of 
    group split differnces by calculating
    pi_between-(pi_within1+pi_within2)/2.
    This assumes that Ne_ancestral is average
    of the two present day Ne_s
    """
    return pwd.subtract(np.diagonal(pwd)/2.,
                          axis=1).subtract(np.diagonal(pwd)/2.,
                                           axis=0)
#    d = np.diagonal(pwd)
#    h = np.repeat([d],len(d),axis=0)
#    v = h.T
#    mx = ( v*(v<h) + h*(v>=h))
#    mn = scipy.stats.gmean([h,v], axis=0)
#    return (pwd - mx)


def get_split_times(pwd, mutation_rate, accessible_size, generation_time=1,
                    groupings=None):
    """
    Uses a pairwise difference matrix to 
    infer group split times. 

    Returns: 
    split_time ... pd.DataFrame of split times
    ne ... dict effective population sizes of each group
    n_samples ... dict of number of haploid samples in each group
    """
    ne = get_ne(pwd, mutation_rate, accessible_size, groupings=groupings)
    
    n_samples = get_samples_per_group(pwd, groupings=groupings, haploid=True)
    
    if groupings is not None:
        pwd = get_group_pwd(pwd, groupings)
    
    
    split_diff = get_split_diff(pwd)
    split_time = split_diff/(2.*mutation_rate*accessible_size/generation_time)
    
    return split_time, ne, n_samples

def get_split_tree(pwd, mutation_rate, accessible_size, generation_time=1,
                    groupings=None, outgroup=None, prune_outgroup=True):
    """
    Uses a pairwise difference matrix to
    infer group split times.
    A neighbour-joining tree is constructed
    from these split times and returned.
    The tree nodes a annotated with branch
    effective population sizes (ne). Ancestral
    ne's are assumed to be average of the child
    branches (but for groupings the average of indivdual
    ne's is taken as basis for this calculation). 

    Returns:
    tree ... pypopgen.modules.treetools.HsTree instance
             (heir of ete3.Tree but with additional functionality,
             e.g., a .plot() method based on matplotlib)
    
    """


    individual_ne = get_ne(pwd, mutation_rate, accessible_size)
    if groupings is not None:
        for k,v in groupings.iteritems():
            if v == outgroup:
                del individual_ne[k]
    else:
        del individual_ne[outgroup]
    if min(individual_ne.values())*2 < max(individual_ne.values()):    
        warnings.warn("Inferred effective population sizes differ by a factor more than 2."
                " The assumptions used to infer split times are not met. " 
                "The tree is likely far off from the truth. Branches with smallest Ne will be far too long. "
                "Here are the estimates: {}".format(str(individual_ne)))
    ne = get_ne(pwd, mutation_rate, accessible_size, groupings=groupings)
    
    n_samples = get_samples_per_group(pwd, groupings=groupings, haploid=True)
    
    if groupings is not None:
        pwd = get_group_pwd(pwd, groupings)
    
    
    split_diff = get_split_diff(pwd)
    split_time = split_diff/(2.*mutation_rate*accessible_size/generation_time)

    #the factor 2 comes from the fact that the distance between two leafes is 2*split_time
    tree = tt.dm_to_tree(2*split_time, outgroup=outgroup, prune_outgroup=prune_outgroup)
    
    
    tree.add_property_to_nodes('ne',ne)
    tree.add_property_to_nodes('n_samples',n_samples)
    
    
    for node in tree.iter_descendants('postorder'):
        if not hasattr(node, 'ne'):
            l,r = node.get_children()
            nes = [l.ne, r.ne]
            for i,n in enumerate([l,r]):
                if n.is_leaf():
                    nes[i] = pwd.loc[n.name, n.name]/(4.*mutation_rate*accessible_size)
                
            node.ne = sum(nes)/2. 
    
    return tree

