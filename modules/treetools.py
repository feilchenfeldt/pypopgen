"""Tools to calculate 4 taxon tests,
such as D-statistic (ABBA-BABA) and 
f-statistic for many samples on a large tree.
"""

# ------------------------------------------------------------------
# Functions to calculate f-statistics (admixture fraction) and jackknife
# ------------------------------------------------------------------


def get_fstat_snpwindow_map_reduce(fn, populations, quadruples,
                                   controlsamples_h3=1, controlsamples_h2=0,
                                   jackknife_window=10000,
                                   chunksize=50000, map_fun=None,
                                   reduce_fun=None):
    def get_fstat_chunk(chunk):
        get_fstat_snpwindow(gen_df, quadruples, populations,
                            controlsamples_h3=controlsamples_h3,
                            controlsamples_h2=controlsamples_h2,
                            jackknife_window=jackknife_window)
