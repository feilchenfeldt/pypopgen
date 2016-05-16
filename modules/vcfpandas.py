"""
Tools that help to load/parse information
from a VCF (variant call format) file into 
pandas data frames and series.
Tools to handle such data.
"""

import re
import subprocess, gzip
import numpy as np
import pandas as pd



def add_info_line(info_dic, line):
    """
    This function parses a vcf header line
    and adds the information in the line
    into info_dic. 
    """
    #split all commas except if in quotes
    r = re.compile('[^,]+".+"[^,]*'
                '|'
                '[^,]+')
    try:
        category, value = line.strip().split('=<')
    except ValueError:
        try: 
            info_dic["uncategorised"].update({line.strip().split('=')[0][2:]:line.strip().split('=')[1]})
        except KeyError:
            try:
                info_dic.update({"uncategorised":{line.strip().split('=')[0][2:]:line.strip().split('=')[1]}})
            except ValueError, e:
                print line.strip().split('=')
                raise e
        return
    category = category[2:]
    tags = r.findall(value[:-1])
    subdic = {k:v for k,v in [t.split('=',1) for t  in tags]}
    try:
        info_dic[category].update({subdic["ID"]:subdic})
    except KeyError:
        info_dic.update({category:{subdic["ID"]:subdic}})
        
def parse_vcf_header(fn):
    """
    This function parses the header of a vcf file
    and returns a list of the table header (last header line)
    and a dictionary of the other header entries.
    """
    info_dic = {}
    fh = gzip.open(fn) if fn[-3:] == ".gz" else open(fn)
    for line in fh:
        if line[:6] == "#CHROM":
            header = line.strip().split('\t')
            header[0] = header[0][1:]
        else:
            add_info_line(info_dic, line)
        if line[0] != '#':
            break
    return header, info_dic

def get_vcf_df(fn, chrom, start=None, end=None, header=None, **read_csv_args):
    """
    This function reads a region from a bgzipped and 
    tabix indexed VCF into a pandas data frame.
    """
    if header is None:
        header, _ = parse_vcf_header(fn)
    region = chrom
    if end is not None and start is None:
        start = 0
    if start is not None:
        region += ':' + str(start)
        if end is not None:
            region += '-' + str(end)

    tabix_stream = subprocess.Popen(['tabix',fn, region], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    vcf_df = pd.read_csv(tabix_stream.stdout, sep="\t", index_col=[0,1], header=False, names = header, **read_csv_args)
    
    return vcf_df

def get_genotype(s):
    """
    Get genotype from vcf genotype string.
    Attention: Only parses biallelic genotypes 
    at the moment.
    """
    gt_str = s.split(':',1)[0]
    if gt_str[:3] in ["0/0","0|0"]:
        return 0
    elif gt_str[:3] in ["1/1","1|1"]:
        return 2
    elif gt_str[:3] in ["0/1","0|1","1|0"]:
        return 1
    else:
        return np.nan
    
def get_info_dic(s):
    """
    Parse the string from the VCF INFO
    column and return it as a dictionary.
    """
    info_tuples = [t.split('=') for t in s.split(';')]
    info_tuples = [t for t in info_tuples if len(t)==2]
    tags = [t for t in info_tuples if len(t)!=2]
#    try:
    d = {k:v for (k,v) in info_tuples}
    d.update({'_tags':tags})
#    except ValueError:
#        d = {}
#        logging.warning("Can not parse info column at {}:{}".format(line[0],line[1]))
    return d


class converters:
    """
    This is just a container for 
    different static converter methods.
    Not meant to be instantiated.
    """
    @staticmethod
    def genotype_converter(samples):
        return {n:get_genotype for n in samples}

def get_genotype_df(fn, chrom, start=None, end=None, header=None):

    if header is None:
        header, _ = parse_vcf_header(fn)
    samples = header[9:]

    vcf_df = get_vcf_df(fn, chrom, start, end, 
                        header=header, 
                        usecols=[0,1] + range(9,len(header)),
                        converters=converters.genotype_converter(samples))
    
    return vcf_df
    
    



