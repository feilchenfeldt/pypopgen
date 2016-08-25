

def get_allele_df_snps(vcf_fn, chrom, samples=None, header=None, start=None,end=None):
    
    def get_alleles_s(ser, sample_ids):
        """
        ser is a series for a given genomic position
        (gen_df.apply(lambda s: get_alleles(s,sample_ids),axis=1)    
        """
        def get_alleles(s):     
            gt_str = s.split(':',1)[0]
            if '|' in gt_str:
                sep = '|' #phased
            else:
                sep = '/' #unphased
            gts = gt_str.split(sep)
            alleles = []
            for ht in gts:
                ht = int(ht)
                if ht == 0:
                    alleles.append(ref)
                else:
                    alleles.append(alts[ht-1])
            if sep == '|' or (len(set(alleles)) == 1):
                return tuple(alleles)
            else:
                return set(alleles)
        
        ref = ser.ix["REF"]
        alts = ser.ix["ALT"].split(',')

        return ser.ix[sample_ids].apply(get_alleles)
    
    if header is None:
        header, _ = vp.parse_vcf_header(vcf_fn)
    if samples is None:
        samples = range(9,len(header))
    raw_allele_df = vp.get_vcf_df(vcf_fn, chrom, start=start, end=end,
                            header=header,usecols=[0,1,3,4]+samples)
    raw_allele_df_snps = raw_allele_df[(raw_allele_df["REF"].apply(len)==1)&\
                                        (raw_allele_df["ALT"].apply(lambda s: all([len(a)==1 for a in s.split(',')]))==1)]
    allele_df = raw_allele_df_snps.apply(lambda s: get_alleles_s(s,samples), axis=1)
    return allele_df

def get_msmc_df(tpl, haplotypes, allele_df, non_callable_s):
    
    def get_msmc_str(allele_s, haplotypes):
        """
        allele_s ... series with allele tuple or set for each indiviual
        haplotypes ... list with, for each indidual 
                        0.. if haplotype 0 should be used
                        1.. if haplotype 1 should be used
                        2.. if both haplotypes should be used
        number of alleles extracted = n_samples + n_samles * sum_i(haploypes[i]==2)
        returns msmc allele string
        """
        assert len(allele_s) == len(haplotypes), "len {} != {}".format(allele_s,haplotypes)
        extracted_alleles_per_individual = []
        for alleles, haplotype in zip(allele_s,haplotypes):
            if haplotype == 2:
                haplotype = [0,1]
            if isinstance(alleles, set):
                #unphased
                allele_permutations = [c for c in itertools.permutations(alleles)]
                used_alleles = ["".join(np.array(t)[haplotype]) for t in allele_permutations]
            else:
                used_alleles = ("".join(np.array(alleles)[haplotype]),)
            extracted_alleles_per_individual.append(used_alleles)
        #for each unphased sample, add all possible phased combinations
        msmc_str = ",".join(["".join(t) for t in itertools.product(*extracted_alleles_per_individual)])
        return msmc_str
    
    def get_n_callable(int_s,non_callable_s):
        filtered_int = non_callable_s.ix[int_s.name[0]].loc[int_s['start']:int_s['end']].reset_index()
        n_filered = (filtered_int['end'] - filtered_int['start']).sum()
        n_callable =  int_s['end'] - int_s['start'] - n_filered
        assert n_callable > 0, ("Something wrong with the .bed "
                                "make sure intervals are non-overlapping, "
                                "0-indexed, right-open and SNPs are spanned by intervals. "
                                "snp_interval: {} \n bed_intervals: {} \n ".format(int_s,filtered_int))
        return n_callable
    
    msmc_s = allele_df[list(tpl)].apply(lambda s: get_msmc_str(s, haplotypes),axis=1)
    msmc_s = msmc_s[msmc_s.apply(lambda s: len(set(s))!=1)]
    msmc_s.name = "alleles"
    
    interval_df = pd.DataFrame({'start': [0] + list(msmc_s.index.droplevel(0))[:-1],
                'end': msmc_s.index.droplevel(0)}, index=msmc_s.index)
    
    n_callable_s = interval_df.apply(lambda int_s: get_n_callable(int_s, non_callable_s),
                                                                                    axis=1)
    n_callable_s.name = 'n_callable'
    msmc_df = pd.concat([n_callable_s, msmc_s],axis=1)
    return msmc_df

def save_msmc_input(vcf_fn, non_callable_bed_fn, chrom, tuples, haplotypes_ls, out_prefix):
    """
    """
    assert len(tuples) == len(haplotypes_ls), "len {} != {}".format(tuples,haplotype_ls)
    
    non_callable_df = pd.read_csv(non_callable_bed_fn,
                sep='\t',header=False,
                names=["chrom","start","end"],index_col=[0])
    
    #convert 0-indexed bed format to 1-indexed VCF coordinates
    non_callable_df[['start','end']] = non_callable_df[['start','end']]+1
    non_callable_s = non_callable_df.set_index('start',append=True)['end']
    
    samples = list(set([i for j in tuples for i in j]))
    allele_df = get_allele_df_snps(vcf_fn, chrom, samples=samples, header=None)
    #return allele_df
    
    for  tpl, haplotypes in zip(tuples,haplotypes_ls):
        try:
            msmc_df = get_msmc_df(tpl, haplotypes, allele_df, non_callable_s)
            msmc_df.to_csv(out_prefix+"_msmc_{}_{}_{}.tsv".format("_".join(tpl),"".join([str(h) for h in haplotypes]),chrom),
                        sep='\t',header=False)
        except AssertionError,e:
            return str(e)


