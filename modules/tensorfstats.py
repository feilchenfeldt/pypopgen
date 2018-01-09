"""
Tensor based implementation of f-statistics to 
simultaneously calculate many comparisons.

Works on VCF files.

For memory efficiency, weighted block-jackknifing is
currently done across whole chromosomes.

F-test was compared to previous implementation (which was checked
against Admixtols). Results were extremely close, e.g.: D +- 0.25%

"""

import logging
import numpy as np
import pandas as pd

#local modules
import vcfpandas as vp
import treetools

logger = logging.getLogger()
logging.basicConfig(
    format='%(levelname)-8s %(asctime)s %(filename)  %(message)s')
#logging.basicConfig(format='%(levelname)-8s %(asctime)s %(funcName)20s()  %(message)s')
logger.setLevel(logging.WARNING)

class Test(object):
    x = 'Test'
    
    def test(self, i):
        return i + self.x
    
    def run_parallel(self, rc):
        rc[:].push({'x': 22})
        rc[:].use_cloudpickle()
        lv = rc.load_balanced_view()
        m = lv.map_async(lambda c: x, range(5))
        #return m.result
        return m.result


class Ftest(object):
    """
    This is the base class for f-tests. Specific
    test classes such as Dtest, F3test, F4test, ...
    derive from it.

    Performs tests for all combinations of h1, h2, h3, (h4).
    ------------------
    Parameters:
    vcf_filename : vcf.gz filename. Should be tabix indexed
                 for random access.
    ind_to_pop : dictionary that maps individuals to populations
                if each individual is its on population this can also
                be a list of individuals
    reduce_dim : If true, remove dimensions with length 1 (not implemented).
    """
    ftype = None

    def __init__(self, vcf_filename,
                ind_to_pop=None,
                result_filebase=None,
                reduce_dim=False):
        

        self.vcf_filename = vcf_filename

        try:
            self.samples = ind_to_pop.keys()
        except AttributeError:
            ind_to_pop = {k:k for k in ind_to_pop}
            self.samples = ind_to_pop.keys()

        self.ind_to_pop = ind_to_pop
        self.result_filebase = result_filebase

    @staticmethod
    def get_hap_df(t0, t1):
        """
        Takes 0/1 DataFrames for two haplotpyes
        per individual and combines them to per
        population allele frequencies.
        """
        t0.columns = pd.MultiIndex.from_arrays(
                [t0.columns, [0] * len(t0.columns)])
        t1.columns = pd.MultiIndex.from_arrays(
                [t1.columns, [1] * len(t1.columns)])
        hap_df  = pd.concat([t0, t1], axis=1).sortlevel(axis=1)
        hap_df = hap_df.dropna(axis=0)
        return hap_df

    @staticmethod
    def get_af(hap_df, ind_to_pop):
        if len(hap_df):
            af = hap_df.groupby(level=0, axis=1).mean()
            af = af.groupby(ind_to_pop, axis=1).mean()
        else:
            af = pd.DataFrame(columns=set(ind_to_pop.values()))
        return af

    @staticmethod
    def get_ac(hap_df, ind_to_pop):
        if len(hap_df):
            ac = hap_df.groupby(level=0, axis=1).sum()
            ac = ac.groupby(ind_to_pop, axis=1).sum()
        else:
            ac = pd.DataFrame(columns=set(ind_to_pop.values()))
        return ac

    @staticmethod
    def get_n(hap_df, ind_to_pop):
        """
        Get the number of haplotypes per population.
        """
        if len(hap_df):
            n = hap_df.groupby(level=0, axis=1).apply(lambda df: df.shape[1])
            
            n = n.groupby(ind_to_pop).sum()
            
        else:
            n = pd.Series(index=set(ind_to_pop.values()))
        return n


    @staticmethod
    def fly_reduce_fun(chunk_res, result=None):
        if result is None:
            return chunk_res
        else:
            return result + chunk_res





    def map(self, chromosomes, start=None, end=None,
             map_fun=map, get_result_fun=lambda r:r, chunksize=50000,
             return_result=True, save_result=False):
        """
        chromosomes : list of chromosome names as in the vcf file
                        to run the analysis on
        stat : Start at this position at each chromosome
        end : End at this position at each chromosomes
            If start and end are not given, the whole chromosome
            is used.
        map_fun : The mapping function that will be used
                  to map calculations to chromosomes. 
                  Could for instance be multiprocessing.map_async
                  or ipython parallel map_async. Default: map 
        get_result_fun : The function that receives the result
                  from the output of the mapping function
                  For map this would simply be the identity (default)
                  or for ipython parallel lv.map_async it is
                  lambda r: r.result() (or lambda r: r.result for older versions) 
        chunksize : Number of lines in the input vcf to read per chunk
                    if jackknife_levels is 'chunk' this also determines
                    the block-jackknife block.
        return_result : Whether to return the result for each chromosme
                        to the mapping function. For very large results
                        it can be more stable and memory efficient to 
                        set return_result=False and save_result=True 
        save_result : whether or not to save the result for each chromosome
                      to disk
        """
        assert return_result or save_result
        if save_result:
            assert self.result_filebase is not None
        result_filebase = self.result_filebase
        self.chromosomes = chromosomes
        self.get_result_fun = get_result_fun
        #calc_stat = lambda chunk1, chunk2: self.calc_stat_static(chunk1, chunk2, ind_to_pop, h1s, h2s, h3s, h4s)        
        #params = {'vcf_filename':self.vcf_filename,
        #         'calc_stat':calc_stat, 'ind_to_pop': self.ind_to_pop, 'h1s':self.h1s,
        #         'h2s':self.h2s, 'h3s': self.h3s, 'h4s': self.h4s,
        #         'samples':self.samples,'fly_reduce_fun':self.fly_reduce_fun,
        #         'chunksize':self.chunksize,'mr_haplo_fun':vp.map_fly_reduce_haplo}    
        vcf_filename = self.vcf_filename
        ind_to_pop = self.ind_to_pop
        samples = self.samples
        fly_reduce_fun = self.fly_reduce_fun
        mr_haplo_fun = vp.map_fly_reduce_haplo

        calc_stat = self.get_calc_stat(*self.calc_params)
        
        def calc_fstat_fun(chrom):
            r = mr_haplo_fun(vcf_filename.format(chrom), calc_stat, 
                                        samples_h0=samples, samples_h1=samples,
                                               chrom=chrom, start=start, end=end, 
                                                fly_reduce_fun=fly_reduce_fun,
                                                               chunksize=chunksize)
            if save_result:
                np.save(result_filebase+'_'+str(chrom), r)
            #return_result = True
            if return_result:
                return r 
        self.map_result =  map_fun(calc_fstat_fun, chromosomes)
        return self.map_result

    def progress(self):
        return self.map_result.progress

    def load_result_chrom(self, chrom):
        r = np.load(self.result_filebase+'_'+str(chrom)+'.npy')
        return r

    def load_result(self):
        res = np.array([self.load_result_chrom(c) for c in self.chromosomes])
        return res

    def get_result(self):
        #try:
        #    se.result(1)
        #except TimeoutError:
        #    logger.INFO('Not finished, status is: {}'.format({i:s.count(i) for i in set(s)}))
        #    return
        try:
            res = np.array(self.get_result_fun(self.map_result))
        except AttributeError:
            res = self.load_result()
        if res[0] is None:
            res = self.load_result()
        stat = self.get_stat(res)
        zscores = self.get_zscores(res, stat)

        stat_df = self.get_stat_df(stat, zscores)
        self.stat_df = stat_df

        return stat_df


class PairwiseDiff(Ftest):
    """
    Parameters:
    hs1s : list of sample names to use as h1
    hs2s : list of sample names to use as h2
    vcf_filename : vcf.gz filename. Should be tabix indexed
                 for random access.
    ind_to_pop : dictionary that maps individuals to populations
                if each individual is its on population this can also
                be a list of individuals
    reduce_dim : If true, remove dimensions with length 1 (not implemented).
    """
    ftype = 'pwd'

    def __init__(self, vcf_filename, ind_to_pop, h1s, h2s, **kwa):
        self.h1s = h1s
        self.h2s = h2s
        self.h3s = None
        self.h4s = None

        self.test_fun = calc.pwd
  

        Ftest.__init__(self, vcf_filename, ind_to_pop, **kwa)

    @staticmethod
    def fly_reduce_fun(chunk_res, result=None):
        """
        Function that reduces on the fly by summing
        """
        if result is None:
            return chunk_res
        else:
            return result + chunk_res 

    @staticmethod
    def calc_stat_static(chunk1, chunk2, ind_to_pop, h1s, h2s, *args):
        hap_df = Ftest.get_hap_df(chunk1, chunk2)
        af = Ftest.get_af(hap_df, ind_to_pop)
        if len(af):
            return calc.pwd(af[h1s], af[h2s])
        else:
            return np.zeros((len(h1s),len(h2s)))

    @staticmethod
    def get_calc_stat(*args):
        def calc_stat(chunk1, chunk2):
            return PairwiseDiff.calc_stat_static(chunk1, chunk2, *args)
        return calc_stat

    @staticmethod
    def jackknife(res, i):
        return np.sum(res[np.arange(len(res))!=i, 0], axis=0)/np.sum(res[np.arange(len(res))!=i, 1], axis=0)
        
    @staticmethod
    def get_stat(res):
        return np.sum(res, axis=0)

    @staticmethod
    def get_zscores(res, pwd):
        return None

    @staticmethod
    def get_stat_df_static(pwd, h1s, h2s):
        return pd.DataFrame(pwd, index=h1s, columns=h2s)

    def get_stat_df(self, stat, zscores):
        return self.get_stat_df_static(stat, self.h1s, self.h2s)



class F3test(Ftest):
    """
    Parameters:
    hs3s : list of sample names to use as h3
          This is the branch called 'C' in
            Patterson et al. 2012 Genetics
    hs1s : list of sample names to use as h1
    hs2s : list of sample names to use as h2
    vcf_filename : vcf.gz filename. Should be tabix indexed
                 for random access.
    ind_to_pop : dictionary that maps individuals to populations
                if each individual is its on population this can also
                be a list of individuals
    reduce_dim : If true, remove dimensions with length 1 (not implemented).
    """
    ftype = 'F3'

    def __init__(self, vcf_filename, ind_to_pop, h3s, h1s, h2s, do_drop_self_comparisons=False, **kwa):
        self.h1s = h1s
        self.h2s = h2s
        self.h3s = h3s
        self.do_drop_self_comparisons = do_drop_self_comparisons
  

        Ftest.__init__(self, vcf_filename, ind_to_pop, **kwa)
        
        self.calc_params = (self.ind_to_pop, self.h1s, self.h2s, self.h3s)

    @staticmethod
    def fly_reduce_fun(chunk_res, result=None):
        """
        Function that reduces on the fly by summing
        """
        if result is None:
            return chunk_res
        else:
            return result + chunk_res 

    @staticmethod
    def calc_stat_static(chunk1, chunk2, ind_to_pop, h1s, h2s, h3s, *args):
        hap_df = Ftest.get_hap_df(chunk1, chunk2)
        af = Ftest.get_af(hap_df, ind_to_pop)
        ac = Ftest.get_ac(hap_df, ind_to_pop)
        n = Ftest.get_n(hap_df, ind_to_pop)
        if len(af):
            return calc.f3(ac[h3s], n[h3s], af[h1s], af[h2s])
        else:
            return np.zeros((len(h3s),len(h1s), len(h2s)))

    @staticmethod
    def get_calc_stat(*args):
        def calc_stat(chunk1, chunk2):
            return F3test.calc_stat_static(chunk1, chunk2, *args)
        return calc_stat

    @staticmethod
    def jackknife(res, i):
        return np.sum(res[np.arange(len(res))!=i, 0], axis=0)/np.sum(res[np.arange(len(res))!=i, 1], axis=0)
        
    @staticmethod
    def get_stat(res):
        return np.sum(res, axis=0)

    @staticmethod
    def get_zscores(res, pwd):
        return None

    @staticmethod
    def drop_self_comparisons_static(stat_df, h1s, h2s, h3s):
        dup_sample_indices = [(h3, h1, h2) for h3 in h3s for h1 in h1s for h2 in h2s if h3==h1 or h3==h2 or h1==h2]
        df = stat_df.drop(dup_sample_indices)
        return df

    def drop_self_comparisons(self):
        df = self.drop_self_comparisons_static(self.stat_df, self.h1s, self.h2s, self.h3s)
        self.stat_df_drop = df
        return self.stat_df_drop

    @staticmethod
    def get_stat_df_static(f3s, stat_name, do_drop_self_comparisons, h1s, h2s, h3s):
        df = pd.DataFrame(f3s.flatten(), 
                index=pd.MultiIndex.from_tuples([(h3,h1,h2) for h3 in h3s for h1 in h1s for h2 in h2s]),
                columns=[stat_name])
        df.index.names = ['h3','h1','h2']
        if do_drop_self_comparisons:
            df = F3test.drop_self_comarisons_static(df, h1s, h2s, h3s)
        return df
    

    def get_stat_df(self, stat, zscores):
        return self.get_stat_df_static(stat, self.ftype, self.do_drop_self_comparisons, self.h1s, self.h2s, self.h3s)



class Dtest(Ftest):
    """
    Parameters:
    hs1s : list of sample names to use as h1
    hs2s : list of sample names to use as h2
    hs3s : list of sample names to use as h3
    hs4s : list of sample names to use as h4
    vcf_filename : vcf.gz filename. Should be tabix indexed
                 for random access.
    ind_to_pop : dictionary that maps individuals to populations
                if each individual is its on population this can also
                be a list of individuals
    reduce_dim : If true, remove dimensions with length 1 (not implemented).


     
    reduce_dim : If true, remove dimensions with length 1 (not implemented).
    jackknife_levels : 'chrom' ... weighted block-jackknife across whole chromosomes
                       'chumk' ... block-jackknife across chunks of chunksize snps
                                   This can be very memory intensive. 
    """
    ftype = 'D'

    def __init__(self, vcf_filename, ind_to_pop, h1s, h2s, h3s, h4s, **kwa):
        self.h1s = h1s
        self.h2s = h2s
        self.h3s = h3s
        self.h4s = h4s

        Ftest.__init__(self, vcf_filename, ind_to_pop, **kwa)

        self.calc_params = (self.ind_to_pop, self.h1s, self.h2s, self.h3s, self.h4s)

    @staticmethod
    def fly_reduce_fun(chunk_res, result=None):
        """
        Function that reduces on the fly by summing
        and also implements a counter of chunks
        so that weighted jackknifing can be performed.
        """
        if result is None:
            return (chunk_res[0], chunk_res[1], 1)
        else:
            return (result[0] + chunk_res[0], result[1] + chunk_res[1], result[2] + 1) 


    @staticmethod
    def calc_stat_static(chunk1, chunk2, ind_to_pop, h1s, h2s, h3s, h4s):
        hap_df = Dtest.get_hap_df(chunk1, chunk2)
        af = Dtest.get_af(hap_df, ind_to_pop)
        if len(af):
            return calc.d(af[h1s], af[h2s], af[h3s], af[h4s])
        else:
            return np.zeros((len(h1s),len(h2s),len(h3s),len(h4s))), np.zeros((len(h1s),len(h2s),len(h3s),len(h4s)))

    @staticmethod
    def get_calc_stat(*args):
        def calc_stat(chunk1, chunk2):
            return Dtest.calc_stat_static(chunk1, chunk2, *args)
        return calc_stat

    @staticmethod
    def jackknife(res, i):
        return np.sum(res[np.arange(len(res))!=i, 0], axis=0)/np.sum(res[np.arange(len(res))!=i, 1], axis=0)
        
    @staticmethod
    def get_stat(res):
        d = np.sum(res[:,0], axis=0)*1./np.sum(res[:,1], axis=0)
        return d

    @staticmethod
    def get_zscores(res, d):
        jackknife_estimates = [Dtest.jackknife(res, i) for i in np.arange(len(res))]
        average = np.average(jackknife_estimates, axis=0, weights=res[:,2]*1./np.sum(res[:,2]))
        variance = np.average(1.*(jackknife_estimates - average)**2, axis=0, weights=res[:,2]*1./np.sum(res[:,2])).astype(float)
        try:
            zscores = d * 1. / ( np.sqrt(variance) * np.sqrt(len(jackknife_estimates)-1) )
        except AttributeError, e:
            print variance.shape
            print np.sqrt(5.)
            print variance.max()
            print variance.min()
            #print np.sqrt(variance)
            return variance
            raise e
        return zscores

    @staticmethod
    def get_stat_df_static(stat, zscores, h1s, h2s, h3s, h4s, stat_name):
        stat_s = pd.Series(stat.flatten(), 
          index=pd.MultiIndex.from_tuples([(h1,h2,h3,h4) for \
                                                           h1 in h1s for h2 in h2s for h3 in h3s for h4 in h4s])) 
        stat_s.name = stat_name

        z_s = pd.Series(zscores.flatten(), 
          index=pd.MultiIndex.from_tuples([(h1,h2,h3,h4) for \
                                                           h1 in h1s for h2 in h2s for h3 in h3s for h4 in h4s])) 
        z_s.name = 'Z'
        stat_df = pd.concat([stat_s, z_s], axis=1)

        stat_df.sort_values(stat_name, ascending=False, inplace=True)

        return stat_df      

    def get_stat_df(self,stat, zscores):
        return self.get_stat_df_static(stat, zscores, self.h1s, self.h2s, self.h3s, self.h4s, self.ftype)

    @staticmethod
    def drop_self_comparisons_static(stat_df, h1s, h2s, h3s, h4s):
        dup_sample_indices = [(h1, h2, h3, h4) for h1 in h1s for h2 in h2s for h3 in h3s for h4 in h4s \
                                                                    if h3==h1 or h3==h2 or h1==h2 or h3==h4 or h1==h4 or h2==h4]
        df = stat_df.drop(dup_sample_indices)
        return df

    def drop_self_comparisons(self):
        df = self.drop_self_comparisons_static(self.stat_df, self.h1s, self.h2s, self.h3s, self.h4s)
        self.stat_df_drop = df
        return self.stat_df_drop

    def get_consistent_with_tree(self, ete_tree):
        """
        Get a data frame with the subset
        of tuples that are consistent with a
        given ete tree.
        Parameters:
        ete_tree : ete3 tree object of all samples.
                   Needs to be rooted and include
                   all outgroups..
        
        """
        self.stat_df_consist = treetools.get_consistent_df(self.stat_df, ete_tree)
        return self.stat_df_consist




class F4ratio(Dtest):
    """
    Parameters:
    hs1s : list of sample names to use as h1
    hs2s : list of sample names to use as h2
    hs3s : list of sample names to use as h3
    hs4s : list of sample names to use as h4
    vcf_filename : vcf.gz filename. Should be tabix indexed
                 for random access.
    ind_to_pop : dictionary that maps individuals to populations
                if each individual is its on population this can also
                be a list of individuals
    reduce_dim : If true, remove dimensions with length 1 (not implemented).


     
    reduce_dim : If true, remove dimensions with length 1 (not implemented).
    jackknife_levels : 'chrom' ... weighted block-jackknife across whole chromosomes
                       'chumk' ... block-jackknife across chunks of chunksize snps
                                   This can be very memory intensive. 
    """
    ftype = 'F4ratio'

    def __init__(self, vcf_filename, ind_to_pop, h1s, h2s, h3s, h4s, **kwa):
        
        Dtest.__init__(self, vcf_filename, ind_to_pop, h1s, h2s, h3s, h4s, **kwa)
        
        pop_to_hap = {pop:[] for pop in set(self.ind_to_pop.values())}


        for s, pop in self.ind_to_pop.iteritems():
            pop_to_hap[pop].append((s, 0))
            pop_to_hap[pop].append((s, 1))
        
        self.pop_to_hap = pop_to_hap
    

        self.calc_params = (self.ind_to_pop, self.pop_to_hap, self.h1s, self.h2s, 
                                                                    self.h3s, self.h4s)


    @staticmethod
    def get_af_hap(hap_df, hap_to_pop):
        if len(hap_df):
            af = hap_df.groupby(hap_to_pop, axis=1).mean()
        else:
            af = pd.DataFrame(columns=set(hap_to_pop.values()))
        return af

    @staticmethod
    def calc_stat_static(chunk1, chunk2, ind_to_pop, pop_to_hap, h1s, h2s, h3s, h4s):
        hap_df = F4ratio.get_hap_df(chunk1, chunk2)
        af = F4ratio.get_af(hap_df, ind_to_pop)
        #do the random subsets for each chunk independently
        
        hap_to_pop_ab = {}
    
        for h3 in h3s:
            samples = pop_to_hap[h3]
            sample_idx = np.arange(len(samples))
            try:
                ixa = np.random.choice(sample_idx, len(samples)/2, replace=False)
            except ValueError, e:
                raise e
            ixb = [i for i in sample_idx if i not in ixa]
            hap_to_pop_ab.update({samples[i]: h3 + '_a' for i in ixa})
            hap_to_pop_ab.update({samples[i]: h3 + '_b' for i in ixb})

        af_sub = F4ratio.get_af_hap(hap_df, hap_to_pop_ab)
        if len(af):
            return calc.f4ratio(af[h1s], af[h2s], af[h3s], af_sub[[h3+'_a' for h3 in h3s]], af_sub[[h3+'_b' for h3 in h3s]], af[h4s])
        else:
            return np.zeros((len(h1s),len(h2s),len(h3s),len(h4s))), np.zeros((len(h1s),len(h2s),len(h3s),len(h4s)))

    @staticmethod
    def get_calc_stat(*args):
        def calc_stat(chunk1, chunk2):
            return F4ratio.calc_stat_static(chunk1, chunk2, *args)
        return calc_stat




class calc:
    """
    This is a container for
    the basic functions that do the
    calculations. 
    Not to be instantiated.

    ATTENTION:
    All the functions that use
    einsum produce nans for any product
    with missing data.
    Missing data should be removed beforehand.
    """

    @staticmethod
    def pwd(af1, af2):
        """
        Calculate pairwise differences (pi and dxy).
        
        Input can be np.ndarray or pd.DataFrame.
        Rows are allele frequencies or haplotypes for variants.
        Columns are individuals or populations.

        Rows containing np.nan should be removed beforehand.

        Result:
        Diagonal entries are pi = 2p(1-p).
        Off-diagonal entries are dxy = pq


        """

        return np.einsum('ij,ik->jk',af1, 1-af2) + np.einsum('ij,ik->jk',1-af1, af2)

    @staticmethod
    def pwd0(af_df):
        """
        Calculate pairwise differences (pi and dxy).
        
        Input can be np.ndarray or pd.DataFrame.
        Rows are allele frequencies or haplotypes for variants.
        Columns are individuals or populations.

        Rows containing np.nan should be removed beforehand.

        Result:
        Diagonal entries are pi = 2p(1-p).
        Off-diagonal entries are dxy = pq


        """
        pw_af = np.einsum('ij,ik->jk',af_df, 1-af_df)
        
        try:
            pw_af = pd.DataFrame(pw_af, index = af_df.columns, columns = af_df.columns)
        except AttributeError:
            pass    

        return pw_af + pw_af.T
 
    @staticmethod
    def f3(ac3, n3, af1, af2):
        """
        Calculate the f3 statistic as defined in
        Patterson et al. Genetics 2012. This is 
        the numerator in Patterson's D statistic.

        This corresponds to the unbiased estimator
        in Patterson et al. 2012 Appendix A.

        af3 corresponds to population 'C' in 
        Patterson et al.

        Input is 3 np.ndarrays or pd.DataFrames
        for each of which rows = SNPs,
        columns = allele frequencies,
        with data for all individuals/populations
        which should
        be tested for h1, h2, h3, respectively.

        Genotypes or haplotypes can be converted
        to allele frequencies as 0, 0.5, 1 and 0, 1,
        respectively.

        Output is a three dimensional np.ndarray,
        j x k x l, where the length of the 3 
        axis is the number of individuals in af2,
        af1, af2, respectively. 
        """
        af3 = ac3*1./n3
        dummy3 = np.ones(af3.shape[1])
        dummy1 = np.ones(af1.shape[1])
        dummy2 = np.ones(af2.shape[1])

        return np.einsum('ij,ij,k,l->jkl',af3, (ac3-1)/(n3-1), dummy1, dummy2)  \
                - np.einsum('ij,ik,l->jkl',af3, af1, dummy2) \
                - np.einsum('ij,k,il->jkl',af3, dummy1, af2) \
                + np.einsum('j,ik,il->jkl',dummy3, af1, af2)
        



    @staticmethod
    def f4(af1, af2, af3, af4):
        """
        Calculate the f4 statistic as defined in
        Patterson et al. Genetics 2012. This is 
        the numerator in Patterson's D statistic.

        Input is 4 np.ndarrays or pd.DataFrames
        for each of which rows = SNPs,
        columns = allele frequencies,
        with data for all individuals/populations
        which should
        be tested for h1, h2, h3, h4, respectively.

        Genotypes or haplotypes can be converted
        to allele frequencies as 0, 0.5, 1 and 0, 1,
        respectively.

        Output is a four dimensional np.ndarray,
        j x k x l x m, where the length of the 4 
        axis is the number of individuals in af1,
        af2, af3, af4, respectively. 
        """
        return (np.einsum('ij,ik,il,im->jklm',af1, 1-af2, 1-af3, af4) \
                - np.einsum('ij,ik,il,im->jklm',1-af1, af2, 1-af3, af4) \
                + np.einsum('ij,ik,il,im->jklm',1-af1, af2, af3, 1-af4) \
                - np.einsum('ij,ik,il,im->jklm',af1, 1-af2, af3, 1-af4))


    @staticmethod
    def f4ratio_denom(af1, af2, af3a, af3b, af4):
        """
        Calculate the f4 statistic as defined in
        Patterson et al. Genetics 2012. This is 
        the numerator in Patterson's D statistic.

        Input is 4 np.ndarrays or pd.DataFrames
        for each of which rows = SNPs,
        columns = allele frequencies,
        with data for all individuals/populations
        which should
        be tested for h1, h2, h3, h4, respectively.

        Genotypes or haplotypes can be converted
        to allele frequencies as 0, 0.5, 1 and 0, 1,
        respectively.

        Output is a four dimensional np.ndarray,
        j x k x l x m, where the length of the 4 
        axis is the number of individuals in af1,
        af2, af3, af4, respectively. 
        """

        dummy = np.ones(af2.shape[1])

        af3a = np.array(af3a)
        af3b = np.array(af3b)

        return (np.einsum('ij,k,il,im->jklm',af1, dummy, (1-af3a) * (1-af3b), af4) \
                - np.einsum('ij,k,il,im->jklm',1-af1, dummy, af3a * (1-af3b), af4) \
                + np.einsum('ij,k,il,im->jklm',1-af1, dummy, af3a * af3b, 1-af4) \
                - np.einsum('ij,k,il,im->jklm',af1, dummy, (1-af3a) * af3b, 1-af4))

    @staticmethod
    def d_denom(af1, af2, af3, af4):
        """
        Calculate numerator of Patterson's
        D statistic. See
        Patterson et al. Genetics 2012.

        Input is 4 np.ndarrays or pd.DataFrames
        for each of which rows = SNPs,
        columns = allele frequencies,
        with data for all individuals/populations
        which should 
        be tested for h1, h2, h3, h4, respectively.

        Genotypes or haplotypes can be converted
        to allele frequencies as 0, 0.5, 1 and 0, 1,
        respectively.

        Output is a four dimensional np.ndarray,
        j x k x l x m, where the length of the 4 
        axis is the number of individuals in af1,
        af2, af3, af4, respectively. 
        """
        return (np.einsum('ij,ik,il,im->jklm',af1, 1-af2, 1-af3, af4) \
                + np.einsum('ij,ik,il,im->jklm',1-af1, af2, 1-af3, af4) \
                + np.einsum('ij,ik,il,im->jklm',1-af1, af2, af3, 1-af4) \
                + np.einsum('ij,ik,il,im->jklm',af1, 1-af2, af3, 1-af4))
   
    @staticmethod
    def d(af1, af2, af3, af4):
        """
        Calculate Patterson's
        D statistic. See
        Patterson et al. Genetics 2012.

        Input is 4 np.ndarrays or pd.DataFrames
        for each of which rows = SNPs,
        columns = allele frequencies,
        with data for all individuals/populations
        which should
        be tested for h1, h2, h3, h4, respectively.

        Genotypes or haplotypes can be converted
        to allele frequencies as 0, 0.5, 1 and 0, 1,
        respectively.

        Output is a four dimensional np.ndarray,
        j x k x l x m, where the length of the 4
        axis is the number of individuals in af1,
        af2, af3, af4, respectively.
        """
        return calc.f4(af1, af2, af3, af4), calc.d_denom(af1, af2, af3, af4)

    @staticmethod
    def f4ratio(af1, af2, af3, af3a, af3b, af4):
        """
        Calculate numerator and denominator of f4 
        admixture ratio. See
        Patterson et al. Genetics 2012.

        Input is 5 np.ndarrays or pd.DataFrames
        for each of which rows = SNPs,
        columns = allele frequencies,
        with data for all individuals/populations
        which should
        be tested for h1, h2, h3a, h3b, h4, respectively.

        af3a and af3b should have the same dimension.
        
        Genotypes or haplotypes can be converted
        to allele frequencies as 0, 0.5, 1 and 0, 1,
        respectively.

        Output is a four dimensional np.ndarray,
        j x k x l x m, where the length of the 4
        axis is the number of individuals in af1,
        af2, af3, af4, respectively.
        """
        return calc.f4(af1, af2, af3, af4), calc.f4ratio_denom(af1, af2, af3a, af3b, af4)


class convert:
    """
    Container object with 
    unctions to manipulate
    data frames etc.

    Not to be instantiated.
    """
    @staticmethod
    def haplo_to_individual(pwd_haplo):
        """
        Function to summarize haplotypes
        of diploids to get individual
        pairwise differences.
        """
        pw = pwd_haplo.groupby(level=0, axis=0).mean().groupby(level=0,axis=1).mean()
        return pw


  









#def get_sample_to_pop(populations):
#    return {vi:k for k,v in populations.iteritems() for vi in v}
