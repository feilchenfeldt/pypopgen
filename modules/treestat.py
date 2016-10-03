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

# local imports
import haplotools as hap

logger = logging.getLogger()
logging.basicConfig(format='%(levelname)-8s %(asctime)s %(filename) '
                    '%(message)s')
logger.setLevel(logging.WARNING)


# -----------------------------------
# -----Get pairwise differences------
# -----------------------------------


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
