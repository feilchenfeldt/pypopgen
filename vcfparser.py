#!/usr/bin/env python
"""
Different functions to parse a  VCF.
See argparse help.

ATTENTION:
If ever adding a reference to the walker as an attribute of the parser,
then the deepcopy in MultiRegionParallelWalker might make problems.

Todo:
Package walkers and parsers in different files.
Auto import all parser files from a given folder.

Long-Term:
Now we parallelise by region. 
If there are very different numbers of entries
 in some regions, the different chunks can have
very different processing time.
Solution:
Directly read the tabix-index and use this information.


#give tempfiles more meaningful names:
outiflename+same random id for all + chunk id



"""
import sys, os, json, uuid, gzip
import logging, argparse, inspect, copy
import numpy as np
import subprocess
import multiprocessing as mp
logger = logging.getLogger()
logging.basicConfig(format='%(levelname)-8s %(asctime)s  %(message)s')
#logging.basicConfig(format='%(levelname)-8s %(asctime)s %(funcName)20s()  %(message)s')
logger.setLevel(logging.DEBUG)
eu = os.path.expanduser
jn = os.path.join
try:
    import pandas as pd
except ImportError:
    logging.warning('Python pandas could not be imported. Several parsers will not work.')

try:
    import tabix
except ImportError:
    logging.warning('pytabix could not be imported. No support for interval or parallel use.')

# assure that all errors get logged
def excepthook(*args):
  logging.getLogger().error('Uncaught exception:', exc_info=args)
sys.excepthook = excepthook



nucleotides = ['A','C','T','G']
no_var = ['N','.','X']





#---------main object-------------------


class Walker(object):
    """
    Generic walker through delimited text file.
    Applies methods of a parser object to each
    line in the file.

    Input:
    in_fh ... file handle or filename of input delimited
              text file. Can be compressed with gzip or bgzip.
    parser... object instance that derives 
              from Parser class
    Examples:

    """
    def __init__(self,in_fh, parser, sep='\t',
                               id='',
                                    skip_multiple_entries=True,
                                    progress_report_interval=50000,**kwa):
        self.id = id
        if not hasattr(in_fh, 'read'):
            logging.info("Input file has no .read method. Assuming it is a filepath string.")
            extension = os.path.splitext(in_fh)[-1]
            if extension in ['.gz','.bgz']:
                logging.info("{}Input file has extension {}. Opening with gzip.".format(self.id,extension))
                in_fh = gzip.open(eu(in_fh))
            else:
                in_fh = open(eu(in_fh))
        else:
            extension = os.path.splitext(in_fh.name)[-1]
            if extension in ['.gz','.bgz']:
                logging.info("{}Input file has extension {}. Opening with gzip.".format(self.id,extension))
                in_fh = gzip.open(eu(in_fh.name))
        self.in_fh = in_fh
        self.sep = sep
        self.parser = parser
        #test this
        #self.parser.walker = self
        self.parser.in_fh = in_fh
        self.skip_multiple_entries = skip_multiple_entries
        self.progress_report_interval = progress_report_interval
        self.finished = False
        self.i = 0
        self.prev_chrom = None
        self.prev_pos = -1
        self.multiple_count = 0
        self.print_warning = True

 
    def _split_line(self,line):
        return line.strip().split(self.sep)

    def _yield_split_line(self,fh):
        while True:
            line = self._split_line(fh.next())
            self._report_progress()
            while self._skip_duplicate_line(line) or self._skip_comment(line):
                line = self._split_line(fh.next())
                self._report_progress()
            yield line



    def _header_line_parser(self,line):
        if self.parser.header_fun is not None:
            self.parser.header_fun(line)


    def _skip_duplicate_line(self,line):
        chrom = line[0]
        pos = int(line[1])
        if chrom == self.prev_chrom:
            assert pos >= self.prev_pos, "vcf positions not in "\
                                            "ascending order at: {}:{},{}".format(chrom,self.prev_pos,pos)
            if pos == self.prev_pos:
                self.multiple_count += 1
                if self.multiple_count > 10 and self.print_warning:
                    logging.warning("Omitting further multiple entry warnings.")
                    self.print_warning = False
                if not self.skip_multiple_entries:
                    if self.print_warning:
                        logging.warning("Multiple entries for pos {}:{}. "
                                  "Keeping all entries.".format(chrom,pos))
                    return False
                else:
                    if self.print_warning:
                        logging.warning("{}Warning, multiple entries for pos {}:{}. "
                              "Skipping all but the first.".format(self.id,chrom,pos))
                    return True
        self.prev_chrom = chrom
        self.prev_pos = pos
        return False


    def _report_progress(self):
        self.i += 1
        if self.i % self.progress_report_interval == 0:
            logging.info("{} Parsed {} lines: {} - {}".format(self.id,self.i,self.prev_chrom,self.prev_pos))

    def _skip_comment(self,line):
        if line[0][0] == "#":
            logging.warning("Skipping comment line in vcf: {}".format(line))
            return True
        else:
            return False


    def parse_header(self):
        logging.info("Parsing header.")
        for line in self.in_fh:
            if line[0] == '#':
                self._header_line_parser(line)
            else:
                break
        for a in self.parser.line_write_attrs:
            try:
                getattr(self.parser,a).flush()
            except AttributeError:
                pass

    def parse(self,fh):
        """
        """
        if self.parser.parse_fun is not None:
            logging.info("{}Parsing vcf body of {}.".format(self.id,fh))
            line_it = self._yield_split_line(fh)
            for d in line_it:
                #logging.debug("hallo")
                self.parser.parse_fun(d)
            #not perfect implementation, prev_pos is not necessarily updated 
            #in children if _skip_duplicate_line is overidden
            logging.info("{}Finished: {} lines at {} {}".format(self.id,self.i,self.prev_chrom,self.prev_pos))

    def cleanup(self):
        """
        """
        #logging.debug("Before cleanup.")
        if self.parser.cleanup_fun is not None:
            logging.info("Starting cleanup.")
            self.parser.cleanup_fun()
        #Not sure why the following was needed,
        #but it causes error when output_fun tries to write to
        #the line write file (maybe just flush instead?)
        #for a in self.parser.line_write_attrs:
        #    try:
        #        getattr(self.parser,a).close()
        #    except AttributeError, e:
        #        logging.warning("Failed to close {}. {}".format(a,str(e)))

    def output(self):
        if self.parser.output_fun is not None:
            logging.info("Creating output.")
            self.result = self.parser.output_fun()
        else:
            self.result = None
        self.finished = True

    def run(self):
        self.parse_header()
        self.parse(self.in_fh)
        self.cleanup()
        self.output()
        logging.info("Run finished.")


class SerialWalker(Walker):
    """
    Walk trough several regions serially.
    Same arguments as Walker, with the addition:
    Input:
    in_fh ... file handle of vcf 
               (must be tabixed and opened with gvcf)
    intervals ... list of string intervals, in Samtools format,
                                            such as Chr1:1-5000
    Note that SerialWalker does not support seperator specification,
    since pytabix does the line splitting automatically.
    """
    def __init__(self,in_fh,parser,intervals,auto_tabix=False, **kwa):
        super(SerialWalker, self).__init__(in_fh, parser, **kwa)
        self.intervals = [self._parse_interval(i) for i in intervals]
        self.auto_tabix = auto_tabix
        self._tabix_init()

    #@staticmethod
    #def str_to_interval(strg):
    #   strg

    def _tabix_init(self):
        self.tabix_fh = tabix.open(self.in_fh.name)
        self._check_tabix()

    def _check_tabix(self):
        try:
            if None in self.intervals[0]:
                self._query_tabix(self.intervals[0][0])
            else:
                self._query_tabix(self.intervals[0])
        except tabix.TabixError, e:
            logging.warning("Tabix raised error: {}".format(str(e)))
            if self.auto_tabix:
                logging.warning("Trying to (re-)index file. This can take a while.")
                p = subprocess.Popen(['tabix','-p','vcf','-f',self.in_fh.name],
                                      stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                out, err = p.communicate()
                if p.returncode != 0:
                    if "was bgzip used to compress this file" in err:
                        base, extension  = os.path.splitext(self.in_fh.name)
                        if extension in ['.gz','.bgz']:
                            logging.warning("File seems not to be compressed with bgzip but ends in .gz or.bgz. "
                                             "Trying to decompress. This can take a while.")

                            p = subprocess.Popen(['gzip','-d','-f',self.in_fh.name],
                                          stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                            out, err = p.communicate()
                            if p.returncode != 0:
                                logging.error(err)
                                raise e
                            name  = base
                            logging.warning("Trying to compress. This can take a while.")
                        else:
                            logging.warning("File seems not to be compressed with bgzip. "
                                             "Trying to compress. This can take a while.")
                            name = self.in_fh.name

                        p = subprocess.Popen(['bgzip','-f',name],
                                      stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                        out, err = p.communicate()
                        if p.returncode != 0:
                            logging.error(err)
                            raise
                        logging.warning("Trying to index file. This can take a while.")
                        self.in_fh = gzip.open(name+'.gz')
                        p = subprocess.Popen(['tabix','-p','vcf','-f',self.in_fh.name],
                                      stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                        out, err = p.communicate()
                        if p.returncode != 0:
                            logging.error(err)
                            raise
                    else:
                        logging.error("Reindexing failed with unhandeled error: {}".format(err))
                        raise e
                self.tabix_fh = tabix.open(self.in_fh.name)
                try:
                    self._query_tabix(intervals[0])
                except tabix.TabixError, e:
                    logging.error("Failed to auto-tabix input file.")
                    logging.error("Is the interval {} in the vcf? ".format(self.intervals[0]))
                    raise e
                logging.info("Auto-tabix successful.")
            else:
                logging.error("Is file compressed with bgzip and does it have tabix index? "
                                " If not, produce it or run with flag/argument auto_tabix.")
                logging.error("Is the interval {} in the vcf? ".format(self.intervals[0]))
                raise e


    def _parse_interval(self,interval):
        try:
            chrompos = interval.split(":")
            chrom = chrompos[0]
            try:
                startend = chrompos[1].split('-')
                start = int(startend[0]) - 1 #numeric intervals in pytabix start at start+1
                try:
                    end = int(startend[1])
                except:
                    end = None
            except IndexError:
                start = None
                end = None
        except AttributeError:
            chrom = interval[0]
            start = interval[1]
            end = interval[2]
        if start is None:
            start = 0
        return [chrom, start, end]


    def _query_tabix(self,interval):
        if interval[1] == 0 and interval[2] is None:
            interval = interval[0]
        try:
            return self.tabix_fh.querys(interval)
        except TypeError:
            try:
                return self.tabix_fh.query(*interval)#this is 1-indexed. the first entry is NOT included!
            except TypeError, e:
                logging.error("Interval: {}: Chromosome must be string and position integer."
                                               " e.g. ('Chr1',1000000,1005000). "
                                "Alternatively use string 'Chr:1000000-1005000'".format(interval))
                raise e
            except tabix.TabixError, e:
                logging.error("Is the interval {} in the vcf? ".format(interval))
                raise e

    def _split_line(self,line):
        """
        Tabix handles the split automatically.
        """
        return line


    def parse(self, fh):
        """
        This is a hacky solution for the fact that
        tabix includes positions outside of the desired interval
        if the positions are deletions and the reference allele overlaps
        with the desired interval.
        """
        if self.parser.parse_fun is not None:
            logging.info("{}Parsing vcf body of {}.".format(self.id,fh))
            line_it = self._yield_split_line(fh)
            logging.debug("Chunk:{}".format(self.parser.chunk))
            for d in line_it:
                #attention: pytabix behaves different from tabix
                #in that the starting index is not in the interval!
                if self.parser.chunk[2] is not None: 
                    if int(d[1]) > self.parser.chunk[1] and int(d[1]) <= self.parser.chunk[2]:
                        self.parser.parse_fun(d)
                    else:
                        logging.warning("{}:Skipping {} cause int(d[1]) >= self.parser.chunk[1] = {}"
                                      "and int(d[1]) < self.parser.chunk[2] = {}".format(
                                          self.parser.chunk,d[1],
                                int(d[1]) >= self.parser.chunk[1],int(d[1]) < self.parser.chunk[2]))
                else:
                    self.parser.parse_fun(d)
            #not perfect implementation, prev_pos is not necessarily updated 
            #in children if _skip_duplicate_line is overidden
            logging.info("{}Finished: {} lines at {} {}".format(self.id,self.i,self.prev_chrom,self.prev_pos))

    def run_no_output(self):
        self.parse_header()
        if self.parser.parse_fun is not None:
            for interval in self.intervals:
                fh = self._query_tabix(interval)
                self.parser.chunk = interval
                self.parser.chrom = interval[0]
                self.parser.pos = interval[1]
                self.parse(fh)
        else:
            logging.warning("No vcf body parse function supplied.")
        self.cleanup()

    def run(self):
        self.run_no_output()
        self.output()
        logging.info("Run finished.")




class ParallelWalker(SerialWalker):
    """
    chunk ...        if chunk is False, intervals are not divided into chunks, i.e.,
                     multiprocessing runs with min(ncpus,len(intervals))
    chunk_factor ... Multiply ncpus by this factor to get the number of chunks.
                     Larger value will imply overhead,
                     but helps to use all processors if some chromosomes are very short.
    """
    def __init__(self, in_fh, parser, intervals=None, ncpus='auto', auto_tabix=False,
                                                            tmp_dir='.',chunk=True,
                                                             chunk_factor=2, **kwa):
        if intervals is None:
            #setup without any interval specific things
            super(SerialWalker, self).__init__(in_fh, parser, **kwa)
            self.intervals = None
            self.auto_tabix = auto_tabix
        else:
            super(ParallelWalker, self).__init__(in_fh, parser,intervals, **kwa)
        self.chunk = chunk
        assert hasattr(self.parser,'reduce_fun'), \
                                            ("Parser {} has no reduce_fun method. "
                                             "No support for parallel execution."\
                                                .format(self.parser.__class__.__name__))
        if self.parser.reduce_fun is None:
            logging.warning("Reduce function is None. Is this really intended?")
        self.kwa = kwa
        if ncpus == 'auto':
            ncpus = mp.cpu_count()
        self.ncpus = ncpus
        self.tmp_dir = os.path.expanduser(tmp_dir)
        if chunk:
            logging.debug("cpus: {}, chunk_factor: {}".format(self.ncpus,chunk_factor))
            self.n_chunks = chunk_factor * self.ncpus
            self.missing = True
            if self.intervals is None:
                self.replace_missing = self.replace_missing_no_interval
                logging.info("No intervals specificed for parallel parsing. "
                             "Trying to infer intervals by searching vcf header. "
                                "Hence, tabix checks will happen at runtime only.")
            else:
                self.intervals = [self._parse_interval(i) for i in self.intervals]
                if not any([i[2] is None for i in self.intervals]):
                    self.missing = False
            if self.missing:
                self._header_line_parser = self._header_line_parser_search_contig_len
        else:
            self.ncpus = min(len(self.intervals),self.ncpus)
            logging.info("Chunking is turned off. {} intervals specified. "
                                                "Parallelising by region. "
                     "Using {} processes.".format(len(self.intervals),self.ncpus))
        self.contic_dic = {}




    def _header_line_parser_search_contig_len(self,line):
        if line[:9] == '##contig=':
                    dic = get_header_line_dic(line)
                    self.contic_dic.update({dic['ID']:int(dic['length'])})
        if self.parser.header_fun is not None:
            self.parser.header_fun(line)

    def replace_missing_no_interval(self):
        logging.info("No intervals given, considering all contigs given in VCF header.")
        self.intervals = [(k,1,v) for k,v in  self.contic_dic.iteritems()]
        super(ParallelWalker, self)._tabix_init()

    def replace_missing(self):
        for interval in self.intervals:
            if interval[2] is None:
                interval[2] = self.contic_dic[interval[0]]
                logging.debug("Inferred interval end from contig info in header: {}".format(interval[2]))
            if interval[1] is None:
                interval[1] = 0
                logging.debug("Assuming interval start to be zero")

    def get_chunks(self):
        logging.debug('Getting junks. Intervals:{} n_chunks: {}'.format(self.intervals,self.n_chunks))
        chunks = self.intervals
        for _ in range(2000):
            lengths = [i[2]-i[1] for i in chunks]
            idx = np.argmax(lengths)
            longest = chunks[idx]
            midpoint = longest[1]+(longest[2]-longest[1])/2
            left_chunk = [longest[0],longest[1],midpoint]
            right_chunk = [longest[0],midpoint,longest[2]]
            chunks[idx] = left_chunk
            chunks.insert(idx+1,right_chunk)
            if len(chunks) >= self.n_chunks:
                break
        else:
            logging.warning("Chunk split algorihm reached hard limit at 2000 splits. Won't create more chunks.")
            logging.warning("Chunks used: {}".format(chunks))
        #chunks.sort()
        logging.debug("Chunks used: {}".format(chunks))
        return chunks

    def setup_subwalkers(self):
        subwalkers = []
        temp_fns = []
        logging.debug("Setting up subwalkers.")
        for i, chunk in enumerate(self.chunks):
            logging.debug("Before copy.")
            #subparser = copy.deepcopy(self.parser)
            subparser = copy.copy(self.parser)
            logging.debug("After copy.")
            subparser.header_fun = None #don't parse header in subparser
            #DON'T do the following, parser must be agnostic to walker, don't pass info to parser
            subparser.chrom = chunk[0] #make chunk information acessible in subparser
            subparser.pos = chunk[1]
            parse_source = inspect.getsource(subparser.parse_fun)
            for a in subparser.line_write_attrs:
                tmp_fn = jn(self.tmp_dir,a+'_'+ '_'.join([str(c) for c in chunk]) +'_'+str(uuid.uuid4()) + ".tmp")
                temp_fns.append(tmp_fn)
                setattr(subparser,a+'_fn',tmp_fn)
            subwalkers.append(SerialWalker(self.in_fh, subparser, [chunk],
                                         id="Chunk {}: ".format(chunk), **self.kwa))
        self.subwalkers = subwalkers
        self.temp_fns = temp_fns

    def reduce(self):
        if self.parser.reduce_fun is not None:
            logging.info("Starting reduce step.")
            self.parser.reduce_fun([p for p in self.subparsers])
        else:
            logging.warning("No reduce function given. Won't do anything.")

    def output(self):
        super(ParallelWalker, self).output()
        for a in self.parser.line_write_attrs:
            try:
                getattr(self.parser,a).close()
            except AttributeError:
                logging.warning("Failed to close {}.".format(a))


    def del_temp_files(self):
        #logging.info("Removing temp files.")
        #while self.temp_fns:
        #    os.remove(self.temp_fns.pop())
        pass

    def run_parser(self, i):
        s = self.subwalkers[i]
        for a in s.parser.line_write_attrs:
            setattr(s.parser,a,open(getattr(s.parser,a+'_fn'),'w'))
        s.run_no_output()
        return s.parser

    def run(self):
        self.parse_header()
        logging.debug('Run function after parsing header.')
        if self.chunk:
            if self.missing:
                self.replace_missing()
            self.chunks = self.get_chunks()
        else:
            self.chunks = self.intervals
        logging.debug('Before setting up subwalkers')
        self.setup_subwalkers()
        subparsers = parmap(self.run_parser,range(len(self.subwalkers)),self.ncpus)
        #this is a hacky solution acconting for the fact that open files cannot be
        #pickled and sent to the child processes by multiprocessing
        for p in subparsers:
            for a in p.line_write_attrs:
                #attention the files are opened in read mode now
                setattr(p,a,open(getattr(p,a+'_fn'),'r'))
        self.subparsers = subparsers
        self.reduce()
        logging.info("Creating output.")
        self.output()
        self.del_temp_files()
        logging.info("Run finished.")



#--------------SUPPORT FUNCTIONS-------------------------

def get_walker(in_fh,parser,intervals=None,ncpus=None,**kwa):
    if ncpus is None or ncpus <= 1:
        if intervals is None:
            return Walker(in_fh,parser,**kwa)
        else:
            if not intervals:
                logging.warning("Intervals given but empty, Walker won't do anything.")
            return SerialWalker(in_fh, parser, intervals,**kwa)
    else:
        return ParallelWalker(in_fh,parser,intervals=intervals,ncpus=ncpus,**kwa)

#parallel support

def fun(f,q_in,q_out):
    while True:
        i,x = q_in.get()
        if i is None:
            break
        q_out.put((i,f(x)))

def parmap(f, X, nprocs):
    q_in   = mp.Queue(1)
    q_out  = mp.Queue()

    proc = [mp.Process(target=fun,args=(f,q_in,q_out)) for _ in range(nprocs)]
    for p in proc:
        p.daemon = True
        p.start()

    sent = [q_in.put((i,x)) for i,x in enumerate(X)]
    [q_in.put((None,None)) for _ in range(nprocs)]
    res = [q_out.get() for _ in range(len(sent))]

    [p.join() for p in proc]

    return [x for i,x in sorted(res)]


#---------support functions used by several parsers------------

def get_header_line_dic(line):
    import re
    try:
        #split all commas except if in quotes
        r = re.compile('[^,]+".+"[^,]*'
                '|'
                '[^,]+')
        tags = r.findall(line.strip().split('=<')[1][:-1])
        dic = {k:v for k,v in [t.split('=') for t  in tags]}
    except ValueError:
        print line.strip().split('=<')[1][:-1].split(',')
        print  [t.split('=') for t  in line.strip().split('=<')[1][:-1].split(',')]
    return dic


def get_012(gt_str):
    if gt_str[:3] in ["0/0","0|0"]:
        return "0"
    elif gt_str[:3] in ["1/1","1|1"]:
        return "2"
    elif gt_str[:3] in ["0/1","0|1","1|0"]:
        return "1"
    elif gt_str[:3] == "./.":
        return "N"
    else:
        raise ValueError("Unsupported genotype " + gt_str)

def revert_gt(gt):
    if gt == '0':
        gt = '2'
    elif gt == '2':
        gt = '0'
    return gt


def get_AA(info_str):
    aa = info_str.split("AA=")[1][0]
    return aa

def add_to_countdic(dic,key):
    try:
        dic[key] += 1
    except KeyError:
        dic[key] = 1

def sum_countdics(countdics):
    tot_countdic = {}
    for countdic in countdics:
        for k,v in countdic.iteritems():
            try:
                tot_countdic[k] += v
            except KeyError:
                tot_countdic[k] = v
    return tot_countdic


def get_info_dic(line):
    info_tuples = [t.split('=') for t in line[7].split(';')]
    info_tuples = [t for t in info_tuples if len(t)==2]
    tags = [t for t in info_tuples if len(t)!=2]
#    try:
    d = {k:v for (k,v) in info_tuples}
    d.update({'_tags':tags})
#    except ValueError:
#        d = {}
#        logging.warning("Can not parse info column at {}:{}".format(line[0],line[1]))
    return d

##--------------------Parsers-----------------------


void_parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
subparsers = void_parser.add_subparsers(dest='parser')


parser_classes = {}
def register(cls):
    parser_classes[cls.__name__] = cls



class MetaParser(type):
    """
    This meta-class handles the creation of subparsers
    for each new parser class on class definition.
    """
    def __new__(cls, clsname, bases, attrs):
        newclass = super(MetaParser, cls).__new__(cls, clsname, bases, attrs)
        if 'register' not in attrs.keys():
            attrs['register'] = True
        if attrs['register']:
            register(newclass)
            new_subparser = subparsers.add_parser(clsname, description=newclass.__doc__)
            args = getattr(newclass,'args')
            newclass.line_write_attrs = []
            newclass.end_write_attrs = []
            newclass.line_read_attrs = []
            if args is not None:
                for arg, pars in args.iteritems():
                    new_subparser.add_argument("--"+arg,**pars)
                    #if argument is a file open in write/read mode 
                    #and it is used in parse_fun --> add it to line_write_attr/ line_read_attrs
                    try:
                        t = pars['type']
                        try:
                            mode = t._mode
                        except AttributeError:
                            try:
                                mode = t.mode
                            except AttributeError:
                                continue
                        if  newclass.parse_fun is not None:
                            parse_code = inspect.getsource(newclass.parse_fun)
                            if arg in parse_code:
                                if mode == 'w':
                                    newclass.line_write_attrs.append(arg)
                                elif mode == 'r':
                                    newclass.line_read_attrs.append(arg)
                                else:
                                    logging.warning("Arg {} of parser class {} has _mode or mode "
                                                    "attribute but mode is not 'w' or 'r'. "
                                                    "This will cause problems in parallel mode "
                                                    "if there is any read or write "
                                                    "happening within the parse method.".format(arg,clsname))
                            #check for files that are written in the end
                            elif mode == 'w' and newclass.output_fun is not None:
                                parse_code = inspect.getsource(newclass.output_fun)
                                if arg in parse_code:
                                    newclass.end_write_attrs.append(arg)
                    except KeyError:
                        pass

            #Add original Parser init to Child init.
            def parserinit(init):
                def newinit(self,**kwa):
                    Parser.__init__(self,**kwa)
                    init(self,**kwa)
                return newinit
            setattr(newclass,'__init__',parserinit(newclass.__init__))
        return newclass


class Parser(object):
    """
    This is the basic parser object.
    Not to be used directly.
    All parsing tools should derive from this class.
    Parsing tools are supplied to a walker
    to be used to parse the lines of the file.
    """
    __metaclass__ = MetaParser
    args = None
    register = False
    def __init__(self,**kwa):
        known_args = self.__class__.args if self.__class__.args is not None else {}
        for arg in kwa:
            assert arg in known_args, ("Unknown argument {},"
                                       "possible args are {}".format(arg,known_args.keys()))
        #this is not really necessary, but useful if line_write_attr are changed in instance
        self.line_write_attrs = copy.copy(self.__class__.line_write_attrs)
        self.end_write_attrs = copy.copy(self.__class__.end_write_attrs)
        self.line_read_attrs = copy.copy(self.__class__.line_read_attrs)
        for arg in known_args:
            try:
                a = kwa[arg]
                if a is None:
                    try:
                        nargs = known_args[arg]['nargs']
                        if nargs in ['*','+']:
                            a = []
                    except KeyError:
                        pass
            except KeyError:
                try:
                    a = known_args[arg]['default']
                    logging.info("Argument {} not supplied, using default {}.".format(arg,a))
                except KeyError:
                    req = False
                    try:
                        if known_args[arg]['required']:
                            req = True
                    except KeyError:
                        pass
                    if not req:
                        a = None
                        try:
                            nargs = known_args[arg]['nargs']
                            if nargs in ['*','+']:
                                a = []
                        except KeyError:
                            pass
                    elif arg in self.end_write_args:
                        a = None
                        logging.info("Argument {} not supplied. It looks like a file that is used for final output. "
                                     "Setting it None and assuming that method output_fun is returning to variable. ")
                    else:
                        raise TypeError("Argument {} not supplied but required.".format(arg))
            setattr(self,arg,a)


        for attr in self.line_write_attrs:
            fh = getattr(self,attr)
            if not hasattr(fh, 'write'):
                logging.info("{} has no .write method. Assuming it is a filepath string: {}.".format(attr,fh))
                fh = open(eu(fh),'w')
            base, extension = os.path.splitext(fh.name)
            if extension in ['.gz','.bgz']:
                logging.info("{} ends in {}, compressing with tabix on the fly.".format(attr, extension))
                base1, ext1 = os.path.splitext(base)
                if ext1 and ext1[1:] in['gff', 'bed', 'sam', 'vcf', 'psltbl', 'gvcf']:
                    logging.info("{} has ending {}, supposing that it is reference ordered data."
                                                                    " Will create tabix index.".format(attr,ext1))
                    index = 'vcf' if ext1[1:] == 'gvcf' else ext1[1:]
                else:
                    logging.info("{} has ending {}. Not recognised."
                                     " Will not create tabix index.".format(attr,ext1))
                    index = False
                setattr(self,attr,TabixWrite(fh,index))
    header_fun = None
    parse_fun = None
    cleanup_fun = None
    output_fun = None

class TabixWrite(object):
    """
    A roughly file-like object,
    That pipes everything that is
    writen to it to bgzip in order to
    compress it.

    Input:
    fh ... file handle
    index ... (bool) Also create tabix index on closing the file.
    """

    def __init__(self,fh, index):
        #idea: inherit all attributes from fh?
        self.fh = fh
        self.name = self.fh.name
        try:
            bgzip = subprocess.Popen(['bgzip','-c'], stdin=subprocess.PIPE, stdout=fh,stderr=subprocess.PIPE)
        except OSError, e:
            logging.error("Error from bgzip output stream, is bgzip installed and in the path?: {}".format(fh))
            raise e
        self.bgzip = bgzip
        self.index = index

    def write(self,string):
        self.bgzip.stdin.write(string)

    def flush(self):
        self.bgzip.stdin.flush()

    def close(self):
        out, err = self.bgzip.communicate()
        self.bgzip.stdin.close()
        if self.bgzip.returncode:
            logging.error("Failed to bgzip {}: {}".format(self.fh.name, err))
            raise IOError
        if self.index:
            logging.info("Creating tabix index for {}".format(self.fh.name))
            tabix = subprocess.Popen(['tabix','-p',self.index,self.fh.name],stderr=subprocess.PIPE)
            out, err = tabix.communicate()
            if self.tabix.returncode or err:
                logging.warning("Failed to create tabix index for {}: {}".format(self.fh.name, err))

class ReduceError(Exception):
    pass

def line_write_reduce_cat(filenames, out_fh):
    command = ["cat"] + filenames
    p = subprocess.Popen(command, stdout=out_fh)
    out, err = p.communicate()
    if p.returncode or err:
        raise ReduceError("Cat reduce step had error and/or non-zero exit status: {}.".format(err))

def line_write_reduce_python(file_handles, out_fh):
    for fh in file_handles:
        for line in fh:
            out_fh.write(line)
        out_fh.flush()

class LineWriteParser(Parser):
    """
    Standard parser that reads from VCF
    and writes to one file line by line 
    (vcf or any other).

    Subclass it and add custom parse method!
    """
    register = False
    args ={
        'out_file':
            {'required':True,
             'type':argparse.FileType('w'),
             'help':"File path to write output to."}
        }
    def reduce_fun(self,selfs):
        try:
            line_write_reduce_cat([s.out_file.name for s in selfs], self.out_file)
            logging.info('Reduced with line_write_reduce_cat on ouf_file.')
        except AttributeError, e:
            try:
                out_stream = self.out_file.bgzip.stdin
                line_write_reduce_cat([s.out_file.name for s in selfs], out_stream)
                logging.info('Reduced with line_write_reduce_cat on bgzip stream.')
            except ReduceError, e:
                logging.warning("Cat reduce step tp bgzip had non-zero exit status: {}."
                            "Trying pytonic reduce.".format(e))
                line_write_reduce_python([s.out_file for s in selfs], self.out_file)
        except Exception, e:
            logging.warning("Cat reduce step failed: {}."
                            "Trying pytonic reduce.".format(e))
            line_write_reduce_python([s.out_file for s in selfs], self.out_file)

class VCFTo012(Parser):
    """
    Extract genotype information into a tsv file.
    Coding 0 for homozygous reference, 1 for heterozygote
    and 2 for homozygous alternative allele.
    Missing genotypes are coded by 'N'
    """
    args ={ 'out_tsv':
                {'required':True,
                'type':argparse.FileType('w'),
                'help':"File path to write output tsv to."}
        }

    def header_fun(self,line):
        if line[1:6] == "CHROM":
            self.out_tsv.write("chrom\tpos\t"+line.split("\t",9)[-1])

    def parse_fun(self,sline):
        gt = map(get_012,sline[9:])
        self.out_tsv.write("\t".join(sline[:2])+"\t"+"\t".join(gt)+"\n")

    def reduce_fun(self,selfs):
        command = ["cat"]+[s.out_tsv.name for s in selfs]
        p = subprocess.Popen(command, stdout=self.out_tsv)
        p.communicate()


class VCFToAncDer012(LineWriteParser):
    args = {'exclude_fixed':{'action':'store_true',
                      'help':'Exclude sites that are fixed '
                           'to outgroup but not seggregating.'},
             'out_file':
                {'required':True,
                'type':argparse.FileType('w'),
                'help':"File path to write output tsv to."}}
    #line write parser uses out_file not out_tsv!!!!!!
    def __init__(self,**kwa):
        self.nts = ['A','C','T','G']

    def header_fun(self,line):
        if line[1:6] == "CHROM":
            self.out_file.write("chrom\tpos\t"+line.split("\t",9)[-1])

    def parse_fun(self, sline):
        pos = int(sline[1])
        chrom = sline[0]
        ref = sline[3]
        alt = sline[4].split(',')
        try:
            aa = get_AA(sline[7])
        except IndexError:
            aa = None
        if sline[6] not in ['.','PASS']: #ignore filterd sites
            pass
        elif aa not in self.nts: #ignore sites with no ancestral state
            pass
        elif len(alt)>1: #ignore multiallelic sites
            pass
        elif len(ref)>1 or len(alt[0])>1: #ignore indels
            pass
        elif aa != ref and \
            (aa != alt[0] and alt[0] in self.nts): #ignore sites where ancestral is 3rd allele
            pass
        elif alt[0] in self.nts or aa != ref: #ingore non-polymporphic non substituted
            gt = map(get_012,sline[9:])
            if aa == alt[0]:
                gt = map(revert_gt, gt)
            self.out_file.write('{}\t{}\t{}\n'.format(chrom,pos,"\t".join(gt)))

class VCFStats(Parser):
    """
    Print some statistics about the variants in the vcf.
    """
    args ={
        'out_fn':
            {'required':True,
             'type':argparse.FileType('w'),
             'help':"File path to write output json to."}
        }


    def __init__(self, **kwa):
        self.var_stats = {}
        self.filters = {}

    def parse_fun(self, line):
        add = lambda s: add_to_countdic(self.var_stats,s)
        addf = lambda s: add_to_countdic(self.filters,s)
        add('total')
        info = get_info_dic(line)
        ref = line[3]
        alt = line[4].split(',')
        pass0 = (line[6] in ['PASS','Pass'])
        if len(alt) > 1:
            add('multiallelic')
        elif len(ref) > 1 or len(alt[0]) > 1:
            add('indels')
        elif alt[0] in nucleotides:
            if float(info['AF'])>0 and float(info['AF'])<1:
                if pass0:
                    add('pass_snp')
                else:
                    add('filter_snp')
                    for f in line[6].split(';'):
                        addf(f)
            else:
                add('snp_non_seg')

        else:
            logging.warning("Unknown variant type at {} - {}: {}".format(line[0],line[1],line[4]))


        try:
            aa = info['AA']
        except KeyError:
            return
        if aa in nucleotides:
            add('ancestral_known')
            if pass0 and len(alt) == 1 and float(info['AF'])>0 and float(info['AF'])<1:
                add('pass_ancestral_known')
                if aa == ref:
                    add('pass_snp_ancestral_is_ref')
                elif aa == alt[0]:
                    add('pass_snp_ancestral_is_alt')
                else:
                    add('pass_snp_ancestral_is_third_allele')




    def reduce_fun(self, selfs):
        self.var_stats = sum_countdics([s.var_stats for s in selfs])
        self.filters = sum_countdics([s.filters for s in selfs])

    def output_fun(self):
        self.var_stats['filters'] = self.filters
        json.dump(self.var_stats, self.out_fn)

class FilterByBed(Parser):
    """
    Add filter tag to sites in intervals of bed file.
    """
    args = {'in_beds':
                      {'required':True,
                       'nargs':'+','type':argparse.FileType('r'),
                       'help':"List of filepathes to the input beds."},
            'filter_names':{'required':True,'nargs':'+',
                            'help':'Filter names corresponding to beds.'},
            'out_vcf':{'required':True,
                       'type':argparse.FileType('w'),
                       'help':'Filepath to output vcf.'}
            }
    def __init__(self,**kwa):
        assert len(self.filter_names)==len(self.in_beds), \
                       "There must be as many filter names as beds."
        self.last_rec = [None for _ in self.in_beds]
        self.seen_chroms = set()

    def header_fun(self,line):
        self.out_vcf.write(line)

    def parse_fun(self,sline):
        def get_rec(fh):
            rec = fh.next().strip().split()
            rec[1] = int(rec[1])
            rec[2] = int(rec[2])
            return rec
        ref = sline[3]
        alt = sline[4].split(',')
        pos = int(sline[1])
        chrom = sline[0]
        self.seen_chroms.update((chrom,))
        for i in range(len(self.last_rec)):
            while self.last_rec[i] is None \
                    or (self.last_rec[i][0]!=chrom and self.last_rec[i][0] in self.seen_chroms) \
                                             or (self.last_rec[i][0]==chrom and self.last_rec[i][2]<pos):
                try:
                    self.last_rec[i] = get_rec(self.in_beds[i])
                except StopIteration:
                    break
            if self.last_rec[i][0] == chrom \
                and self.last_rec[i][1] < pos \
                and self.last_rec[i][2] + 1 > pos:
                if sline[6] in ['.','PASS']:
                    sline[6] = self.filter_names[i]
                else:
                    sline[6] = sline[6] + ','  + self.filter_names[i]
        self.out_vcf.write("\t".join(sline)+'\n')


class GetFilterStats(Parser):
    """
    Count occurences of all combinations
    of filters in the filter column.
    """
    args = {'out_fn':{'required':True,
                      'type':argparse.FileType('w'),
                      'help':"Filename to write to."}}

    def __init__(self,**kwa):
        try:
            pd
        except NameError, e:
            logging.error("{} requires python pandas, but it seems not to be imported.".format(self.__class__.__name__))
            raise e
        self.count_dic = {}

    def parse_fun(self,sline):
        filters = sline[6].split(';')
        filters.sort()
        filters = tuple(filters)
        add_to_countdic(self.count_dic, 'n_sites')
        add_to_countdic(self.count_dic, filters)

    def cleanup_fun(self):
        filter_info = pd.Series(self.count_dic.values(),
                                    index=self.count_dic.keys())
        filter_info.sort(inplace=True, ascending=False)
        self.filter_info = filter_info

    def reduce_fun(self,selfs):
        fi = selfs[0].filter_info
        for ad in selfs[1:]:
            fi = fi.add(ad.filter_info, fill_value=0)
        self.filter_info = fi

    def output_fun(self):
        if self.out_fn is not None:
            self.filter_info.to_csv(self.out_fn,sep='\t')
        else:
            return self.filter_info

class AccessibleGenomeStats(Parser):
    """
    Parse a whole genome (all sites) VCF 
    to get statistics about the number of called
    SNPs and nonsnps. This gives information on the number
    of accessible (i.e., non Filtered, non missing genotype)
    sites both global and per individual.
    """
    args = {"out_filter_count":{'required':True,
                      'type':argparse.FileType('w'),
                      'help':"File path to write count of filtered sites as json to."},
            "out_N_count":{'required':True,
                      'type':argparse.FileType('w'),
                      'help':"Tsv path to write per individual missing genotype count to."},
            "out_N_corr":{'required':True,
                      'type':argparse.FileType('w'),
                      'help':"Tsv path to write cross individual correlations "
                                                "of missing genotypes to."}}

    def __init__(self,**kwa):
        self.sites_dic = {"total":0}#"pass_nosnp":0,
                                #"pass_snp":0,"filter_nosnp":0,"filter_snp":0}


    def header_fun(self,line):
        if line[:6] == '#CHROM':
            self.samples = line.strip().split('\t')[9:]
            self.N_df = pd.DataFrame(0,
                                            columns=self.sites_dic.keys(),
                                            index=self.samples)
            self.Nxy = np.zeros((len(self.samples),len(self.samples)))

    def parse_fun(self, line):
        add = lambda s: add_to_countdic(self.sites_dic,s)

        #sites_dic = {"total":0,"pass_nosnp":0,"pass_snp":0,"filter_nosnp":0,"filter_snp":0}

        add("total")
        ns = np.array([1 if '.' in s.split(':')[0] else 0 for s in line[9:]]) #vector of Ns
        self.N_df["total"] +=  ns
        self.Nxy += np.outer(ns,ns)
        ref = line[3]
        alt = line[4].split(',')
        pass0 = (line[6] in ['PASS','Pass'])

        if pass0:
            category = 'pass_'
        else:
            category = 'filter_'
        try:
            af = float(get_info_dic(line)['AF'].split(',')[0])
            if len(alt) == 1 and len(ref) == 1 and len(alt[0]) == 1 and alt[0] in nucleotides and af>0 and af<1:
                category += 'snp'
            else:
                category += 'nosnp'
        except KeyError:
            category += 'nosnp'

        #test where the discrepancy to SNPstats comes from....
#        try:
#            af = float(get_info_dic(line)['AF'].split(',')[0])
#            if len(alt) != 1: #and len(ref) == 1 and len(alt[0]) == 1 and alt[0] in nucleotides and af>0 and af<1:
#                category += 'altgt1'
#            elif len(ref) != 1:
#                category += 'refgt1'
#            elif len(alt[0])>1:
#                category += 'alt0gt1'
#            elif alt[0] not in nucleotides:
#                category += '?allele'
#            elif af==0:
#                category += 'af0'
#            elif af==1:
#                category += 'af1'
#            else:
#                category += 'snp'
#        except KeyError:
#            category += 'noaf'

        add(category)
        try:
            self.N_df[category] += ns
        except KeyError:
             self.N_df[category] = ns

    def reduce_fun(self, selfs):
        #self.N_df = sum([s.N_df for s in selfs])
        self.N_df = reduce(lambda df0, df1: df0.add(df1, fill_value=0), [s.N_df for s in selfs])
        self.N_df = self.N_df.astype(int)
        self.sites_dic = sum_countdics([s.sites_dic for s in selfs])
        self.Nxy = sum([s.Nxy for s in selfs])

    def output_fun(self):
        N_df = self.N_df
        sites_dic = self.sites_dic
        Nxy = pd.DataFrame(self.Nxy,index=self.samples,columns=self.samples)
        corr = (Nxy-1./sites_dic["total"]*np.outer(N_df["total"],N_df["total"]))/\
                np.sqrt(np.outer(N_df["total"]*(1-1./sites_dic["total"]*N_df["total"]),N_df["total"]*(1-1./sites_dic["total"]*N_df["total"])))
        #try:
        json.dump(sites_dic,self.out_filter_count)
        N_df.to_csv(self.out_N_count,sep='\t')
        corr.to_csv(self.out_N_corr,sep='\t')
        #except KeyError:
        #    logging.info("At least one output filename not found. Returning results.")
        #    return (sites_dic, N_df, corr)



class AddFilterInfo(Parser):
    """
    Adds info header entries for each pair
    in expressions/descriptions iterables to
    in_vcf and writes it to out_vcf
    """
    args = {
        'expressions': {'required':True,
                        'nargs':'+',
                        'help':"List of filter expressions to add to the vcf header."},
        'descriptions': {'required':True,
                        'nargs':'+',
                        'help':'List of filter descriptions to add to the vcf header.'},
         'out_vcf': {'required':True,
                     'type': argparse.FileType('w'),
                      'help':"Filepath to write output vcf to."}}
    def __init__(self,**kwa):
        self.line_counter = 0
        self.filter_found = False
        self.filters_added = False

    def header_fun(self,line):
        self.line_counter += 1
        b = len("##FILTER=<ID=")
        if line[:b] == "##FILTER=<ID=":
            logging.debug("Exisitng filter found: {}".format(line))
            self.filter_found = True
            for j,e in enumerate(self.expressions):
                if line[b:b+len(e)] == e:
                    line = line[:b+len(e)+len(',Description="')] + self.descriptions[j] + '">\n'
                    logging.info("Updating filter expression {}.".format(e))
                    del self.expressions[j]
                    del self.descriptions[j]
                    break
        elif not self.filters_added and (self.filter_found or (line[:6] == '#CHROM')):
            for e,d in zip(self.expressions,self.descriptions):
                self.out_vcf.write('##FILTER=<ID=' + e + ',Description="' + d + '">\n')
                logging.info("Adding filter expression {}.".format(e))
            self.filters_added = True
        self.out_vcf.write(line)



    def cleanup_fun(self):
        self.out_vcf.flush()
        logging.info("Starting to cat vcf body.")
        command = ["tail","-n +" + str(self.line_counter+1),self.in_fh.name]
        p = subprocess.Popen(" ".join(command),shell=True,stdout=self.out_vcf)
        #p.wait()
        p.communicate()

class FiltersToFasta(LineWriteParser):#(LineWriteParser):
    """
    Write a fasta where all filtered sites are set to N.
    """
    args = {'type':{'default':'ref','choices':['ref','alt','anc'],
                     'help':"Specify which base to write to the fasta."
                            "anc required AA tag in INFO column."},
            'line_length':{'type':int,'default':80,
                           'help':'Characters per line in output fasta.'},
            'out_file':
                {'required':True,
                'type':argparse.FileType('w'),
                'help':"File path to write output fasta to."}}
    def __init__(self,**kwa):
        self.chrom = None
        self.pos = None
        if self.type == 'anc':
            self.ancestral_tag = False

    def header_fun(self, line):
        if self.type == 'anc':
            if line[:len('##INFO')] == '##INFO':
                d = get_header_line_dic(line)
                if d['ID']=='AA':
                    self.ancestral_tag = True
            if line[:len("#CHROM")] == '#CHROM':
                assert self.ancestral_tag, ("VCF has no INFO tag 'AA' defined in header, "
                                            "won't run with type 'anc'. ")


    def parse_fun(self, sline):
        chrom = sline[0]
        pos = int(sline[1])
        ref = sline[3]
        filter = sline[6]
        pass0 = filter in ['PASS','Pass']
        #attention, there is no check whether 1 is missing
        #but this is tricky in parallel mode
        #logging.info(self.chunk)
        if self.pos is None:
            self.pos = self.chunk[1]
        if self.pos == 0:
            #print 'writing chrom'
            self.out_file.write('>' + chrom)
        while self.pos + 1 < pos:
            if not self.pos % self.line_length:
                self.out_file.write('\n')
            self.out_file.write('N')
            self.pos += 1
        if self.pos >= pos:
            return
        self.chrom = chrom
        self.pos = pos
        if not pass0:
            allele = 'N'
        elif self.type == 'anc':
            try:
                allele = get_AA(sline[7])
            except IndexError:
                allele = 'N'
        elif self.type == 'ref':
            if len(ref)>1:
                print pos
                ref = ref[0]
            allele = ref
        elif self.type == 'alt':
            allele = alt
        else:
            raise NotImplementedError("Mode {} not implemented.".format(self.type))

        if not (self.pos - 1) % self.line_length:
            self.out_file.write('\n')
        self.out_file.write(allele)
        #self.out_file.write(str(pos)[-1])

    def cleanup_fun(self):
        while self.pos < self.chunk[2]:
            if not (self.pos - 1) % self.line_length:
                self.out_file.write('\n')
            self.out_file.write('N')
            self.pos += 1
    def output_fun(self):
        self.out_file.write('\n')

class FiltersToBed(LineWriteParser):
    """
    Write all filter entries into bed file.
    Only works with all-sites VCF.
    """
    args = {#'one_file':{'type':bool,'action':'store_true','default':'False',
            #         'help':"Write entries for all filter columns in a single file. "
            #                "Otherwise one file per filter"},
            'out_file':
                {'required':True,
                'type':argparse.FileType('w'),
                'help':"File path to write output fasta to."}}
    def __init__(self,**kwa):
        self.filter_open = {}

    def header_fun(self, line):
        if line[:len('##FILTER')] == '##FILTER':
           filter_line_dic = get_header_line_dic(line)
           self.filter_open.update({filter_line_dic['ID']:[None,None,None]})

    def parse_fun(self, sline):
        chrom = sline[0]
        pos = int(sline[1])
        filter = sline[6]
        if filter not in ['.','PASS','Pass']:
            filters = filter.split(";")
            for f, (chrom0, start, end) in self.filter_open.iteritems():
                if chrom0 is not None \
                    and chrom0 != chrom:
                        self.out_file.write("{}\t{}\t{}\t{}\n".format(chrom0, start, end, f))
                        self.filter_open[f] = chrom0, start, end  = [None, None, None]
                if chrom0 == None:
                    if f in filters:
                        self.filter_open[f] = [chrom, pos - 1, pos]
                else:
                    if f in filters:
                        self.filter_open[f][2] = pos
                    else:
                        self.out_file.write("{}\t{}\t{}\t{}\n".format(chrom, start, end, f))
                        self.filter_open[f] = [None, None, None]
        else:
            for f, (chrom0, start, end) in self.filter_open.iteritems():
                if chrom0 is not None:
                    self.out_file.write("{}\t{}\t{}\t{}\n".format(chrom0, start ,end, f))
                    self.filter_open[f] = [None, None, None]



class FiltersToFasta2(LineWriteParser):#(LineWriteParser):
    """
    Write a fasta where all filtered sites are set to N.
    Required whole genome-vcf.
    This function is tricky, because:
        -- sometimes whole genome vcfs start at position 2
        -- how to retain information on new chromosome
            start, when chunking hits precisely a chromosome
            boundary?
    ATTENTION: PARALLEL mode not supported!
    """
    args = {'type':{'default':'ref','choices':['ref','alt','anc'],
                     'help':"Specify which base to write to the fasta."
                            "anc required AA tag in INFO column."},
            'line_length':{'type':int,'default':80,
                           'help':'Characters per line in output fasta.'},
            'out_file':
                {'required':True,
                'type':argparse.FileType('w'),
                'help':"File path to write output fasta to."}}
    parallel = False
    def __init__(self,**kwa):
        self.chrom = ''
        self.pos = 1
        if self.type == 'anc':
            self.ancestral_tag = False

    def header_fun(self, line):
        if self.type == 'anc':
            if line[:len('##INFO')] == '##INFO':
                d = get_header_line_dic(line)
                if d['ID']=='AA':
                    self.ancestral_tag = True
            if line[:len("#CHROM")] == '#CHROM':
                assert self.ancestral_tag, ("VCF has no INFO tag 'AA' defined in header, "
                                            "won't run with type 'anc'. ")


    def parse_fun(self, sline):
        chrom = sline[0]
        pos = int(sline[1])
        ref = sline[3]
        filter = sline[6]
        pass0 = filter in ['PASS','Pass']
        if self.pos == 1 or chrom != self.chrom:
            #print 'writing chrom'
            self.out_file.write('>' + chrom)
            if chrom != self.chrom:
                self.pos = 1
        while self.pos < pos:
            if self.pos % self.line_length == 0:
                self.out_file.write('\n')
            self.out_file.write('N')
            self.pos += 1
        if self.pos >= pos:
            raise UserException()
        self.chrom = chrom
        self.pos = pos
        if not pass0:
            allele = 'N'
        elif self.type == 'anc':
            try:
                allele = get_AA(sline[7])
            except IndexError:
                allele = 'N'
        elif self.type == 'ref':
            if len(ref)>1:
                ref = ref[0]
            allele = ref
        elif self.type == 'alt':
            allele = alt
        else:
            raise NotImplementedError("Mode {} not implemented.".format(self.type))

        self.out_file.write(allele)
        if self.pos % self.line_length == 0:
            self.out_file.write('\n')
        #self.out_file.write(str(pos)[-1])

    def cleanup_fun(self):
        while self.pos < self.chunk[2]:
            if not (self.pos - 1) % self.line_length:
                self.out_file.write('\n')
            self.out_file.write('N')
            self.pos += 1
    def output_fun(self):
        self.out_file.write('\n')


class RemoveLowQualNonVariants(LineWriteParser):
    #args ={
    #    'out_file':
    #        {'required':True,
    #         'type':argparse.FileType('w'),
    #         'help':"File path to write output to."}
    #    }
    def header_fun(self,line):
        self.out_file.write(line)
    def parse_fun(self, sline):
        #generalise!!!!!!
        try:
            if sline[4] == '.' and float(sline[5])>10:
                filters = sline[6].split(';')
                try:
                    filters.remove('LowQual')
                except ValueError:
                    pass
                if not filters:
                    sline[6] = 'PASS'
                else:
                    sline[6] = ';'.join(filters)
        except ValueError:
            pass
            #print sline
        self.out_file.write("\t".join(sline)+'\n')

    #def reduce_fun(self,selfs):
        #try:
        #line_write_reduce_cat([s.out_file.name for s in selfs], self.out_file)
        #except AttributeError:
            #logging.warning("Cat reduce step had non-zero exit status: {}."
            #                "Trying pytonic reduce.".format(err))
    #    line_write_reduce_python([s.out_file for s in selfs], self.out_file)

class AddAncestralFasta(LineWriteParser):
    """
    Add ancestral state from fasta to vcf.
    E.g.,
    Macaque state is taken from vervet-macaque
    alignment and added in the info column with
    tag AA.
    """
    args = {'ancestral_source':{'required':True,
                                'help':"Name of source of ancestral allele info."},
              'ancestral_fasta':{'required':True,
                                 'type':argparse.FileType('r'),
                                 'help':'Filepath of fasta with ancestral state.'},
              'out_file':{'required':True,
                                 'type':argparse.FileType('w'),
                           'help':"Filepath to write output vcf to."}}
    
    def __init__(self,**kwa):
        from pyfasta import Fasta
        self.fasta = Fasta(self.ancestral_fasta.name)
        self.info_parsed = False
        self.info_written = False

    def header_fun(self, line):
        #print self.parse_fun
        if line[:7] == '##INFO=':
             self.info_parsed = True
        elif self.info_parsed and not self.info_written:
            self.out_file.write(\
                                    '##INFO=<ID=AA,Number=1,Type=Character'
                                    ',Description="Ancestral allele as'
                                    ' derived from {}">\n'.format(self.ancestral_source))
            self.info_written = True
        self.out_file.write(line)
        


    def parse_fun(self, sline):
        #print "Hallo"
        aa = self.fasta[sline[0]][int(sline[1])-1]
            #aa = 'X'
        if aa not in ['N','n']:
           sline[7] = sline[7] + ';AA=' + aa
        #logging.info("Writing line {} to {}:".format(sline,self.out_file))
        self.out_file.write("\t".join(sline)+'\n')


#    def reduce_fun(self,selfs):
#        try:
#            line_write_reduce_cat([s.out_file.name for s in selfs], self.out_file)
#            logging.info('Reduced with line_write_reduce_cat on ouf_file.')
#        except AttributeError, e:
#            try:
#                out_stream = self.out_file.bgzip.stdin
#                line_write_reduce_cat([s.out_file.name for s in selfs], out_stream)
#                logging.info('Reduced with line_write_reduce_cat on bgzip stream.')
#            except ReduceError, e:
#                logging.warning("Cat reduce step tp bgzip had non-zero exit status: {}."
#                            "Trying pytonic reduce.".format(e))
#                line_write_reduce_python([s.out_file for s in selfs], self.out_file)
#        except Exception, e:
#            logging.warning("Cat reduce step failed: {}."
#                            "Trying pytonic reduce.".format(e))
#            line_write_reduce_python([s.out_file for s in selfs], self.out_file)


#class FilterGenotypes(LineWriteParser):
#    """
#    Add ancestral state from fasta to vcf.
#    E.g.,
#    Macaque state is taken from vervet-macaque
#    alignment and added in the info column with
#    tag AA.
#    """
#    args = {'ancestral_source':{'required':True,
#                                'help':"Name of source of ancestral allele info."},
#              'ancestral_fasta':{'required':True,
#                                 'type':argparse.FileType('r'),
#                                 'help':'Filepath of fasta with ancestral state.'},
#              'out_file':{'required':True,
#                                 'type':argparse.FileType('w'),
#                           'help':"Filepath to write output vcf to."}}
#....
#
#    def header_fun(self, line):
#        self.out_file.write(line)
#
#
#    def parse_fun(self, sline):
#        for i,gt in enumerate(sline[9:]):
#            if gt[:3] == '0/1':
#                sline[9+i] = './.'+ sline[9+i][3:]
#        self.out_file.write("\t".join(sline)+'\n')




class SNPEFFParser(Parser):
    """
    Extract genotype information into a tsv file.
    Coding 0 for homozygous reference, 1 for heterozygote
    and 2 for homozygous alternative allele.
    Missing genotypes are coded by 'N'
    """
    args ={
        'out_tsv':
            {'required':True,
             'type':argparse.FileType('w'),
             'help':"File path to write output tsv to."}
        }


    CHR_IX =0
    POS_IX =1
    REF_IX =3
    ALT_IX =4
    INFO_IX =7
    AT_START_IX = 9
    AT_LENGTH = 1135
    AT_END_IX = AT_START_IX+ AT_LENGTH
    LYR_START_IX = 1171
    LYR_END_IX = LYR_START_IX +2
    REF = "0"
    ALT = "1"
    DOUBLE_REF = "00"
    DOUBLE_ALT = "11"

    def header_fun(self,line):
        if line[1:6] == "CHROM":
            self.out_tsv.write("chrom\tpos\tref\talt\tlyr\tgenotype\tinfo\n")

    def parse_fun(self,sline):
        output = self._filter_(sline)
        if output is not None:
            self.out_tsv.write("\t".join(output)+"\n")

    def reduce_fun(self,selfs):
        command = ["cat"]+[s.out_tsv.name for s in selfs]
        p = subprocess.Popen(command, stdout=self.out_tsv)
        p.communicate()


    def _filter_(self,fields):
        if fields[self.ALT_IX] == ".":
            return None
        if len(fields[self.REF_IX]) > 1 :
            return None
        lookupMap = self._create_set_(fields[self.AT_START_IX:self.AT_END_IX])
        if len(lookupMap) != 2:
            return None
        alt = fields[self.ALT_IX]
        genotype = "1"
        if len(alt) > 1:
            alt_split = alt.split(",")
            if len(alt_split)  == 1:
                return None
            keys = set(lookupMap.keys())
            keys = sorted(keys)
            genotype = keys[-1][0]
            ix = int(genotype)
            if len(alt_split) >= ix:
                alt = alt_split[ix-1]
                if len(alt) != 1:
                    return None
            else:
                return None
        lookupMap = self._create_set_(fields[self.LYR_START_IX:self.LYR_END_IX])
        lyr = self._get_lyr_allele_(lookupMap)
        return [fields[self.CHR_IX],fields[self.POS_IX],fields[self.REF_IX],alt,lyr,genotype,fields[self.INFO_IX]]



    def _create_set_(self,fields):
        lookupMap = {}
        for item in fields:
            sep = item[1]
            isValid = False
            if sep == '|':
                isValid = True
            elif sep == '/':
                if item[0] != '.' and item[2] != '.':
                    isValid = True
            if not isValid:
                continue
            genotype = item[0]+item[2]
            lookupMap[genotype] = True
        return lookupMap

    def _get_lyr_allele_(self,lookupMap):
        size = len(lookupMap)
        if size == 2:
            return 'S'
        elif size == 1:
            if self.DOUBLE_REF in lookupMap:
                return '0'
            if self.DOUBLE_ALT in lookupMap:
                return '1'
            if self.REF+self.ALT in lookupMap:
                return 'S'
            if self.ALT+self.REF in lookupMap:
                return 'S'
        return 'NA'







def get_argparser():
    argparser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                        description="Parse a Variant Call Format (VCF) file.",
                                                                                 add_help=False)


    argparser.add_argument("--variant",'-V', type = argparse.FileType('r'),
                                #default = '-',
                                 help="Input vcf filepath.")

    argparser.add_argument("--parser",'-P', choices = subparsers.choices.keys(),
                                                                help="Subparser to be used.")

    argparser.add_argument("--intervals",'-L', nargs='*', dest='intervals', action='append',
                            help='Specify intervals to consider e.g. Chr1:1-50000. '
                                 'Input vcf must be compressed with bgzip and indexed '
                                                                         'with tabix.')
    argparser.add_argument("--auto_tabix",action='store_true',
                                help="Automatically compress and/or index on the fly.")
    argparser.add_argument("--ncpus", '-nct',
                            type=int, default=1,
                                  help='Number of processes for parallel parsing. '
                                       'If no intervals are specified with '
                                       'with -L (e.g. -L Chr1:1-5000000), '
                                       'or if intervals lack the end position (e.g. -L Chr1), '
                                       'then we try to infer all regions present in the VCF '
                                       'from the VCF header (tag contig=...). '
                                       'Make the header contig tags are consisten with the '
                                       'body in this case. '
                                       'Parallel parsing is not implemented for all '
                                       'parsers.')
    argparser.add_argument("--no_chunk", action='store_true',
                            help="If no_chunk is set, intervals are not divided into chunks, i.e., "
                                 "multiprocessing runs with min(ncpus,len(intervals))")
    argparser.add_argument("--chunk_factor", default=4,type=int,
                            help="Multiply ncpus by this factor to get the number of chunks. "
                                  "Larger value will imply overhead, "
                                  "but helps to use all processors for the whole runtime "
                                  "when chromosomes or intervals are of unequal length.")
    argparser.add_argument("--temp_dir",
                                 default='.',
                                  help= "Path to folder to write temporary files to. "
                                        "Only used when parallel processing.")
    argparser.add_argument("--skip_multiple_entries",action='store_true',
                         help='Skip all but the first entry for the same site in the VCF.')
    argparser.add_argument('--logging_level','-l',
                        choices=['DEBUG','INFO','WARNING','ERROR','CRITICAL'],
                        default='INFO',
                                                help='Minimun level of logging.')
    argparser.add_argument('--progress_report_interval',
                        type=int, default=400000,
                        help='Number of lines after which to report progress. '
                                            'Only output if logging level >= INFO')
    argparser.add_argument("--help",'-h', action='store_true',
                                                     help="Print help message and exit.")


    return argparser

def parse_args(argparser):
    args, unknown = argparser.parse_known_args()

    if args.help:
        help_str = "\n"+argparser.format_help()
        logger.setLevel(logging.INFO)
    try:
        subargparser = subparsers.choices[args.parser]
        if args.help:
            help_str += "\n\n-----------Help for parser {} ------------\n\n".format(args.parser) \
                                                                         + subargparser.format_help()
    except KeyError, e:
        if args.help:
            pass
        elif args.parser is None:
            argparser.error("argument --parser/-P is required")
        else:
            raise e

    if args.help:
       logging.info(help_str)
       sys.exit(0)

    if args.variant is None:
        argparser.error("argument --variant/-V is required")

    #if args.variant is sys.stdin:
    #    logging.warning("argument --variant/-V not supplied, trying to read from stdin")

    sub_args = subargparser.parse_args(unknown)
    return args, sub_args

def parse(args,sub_args):
    import select
    import time
    logger.setLevel(getattr(logging,args.logging_level))
    if args.intervals is not None:
        args.intervals = [a for b in args.intervals for a in b]

    try:
        assert select.select([args.variant,],[],[],0.0)[0], "Input vcf has no data."
    except TypeError: #the above check does not work for tabix
        pass

    parser_class = parser_classes[args.parser]
#    logging.warning("INTERVAL INPUT: {}".format(args.intervals))
    parser = parser_class(**{arg:getattr(sub_args,arg) for arg in vars(sub_args)})


    walker = get_walker(args.variant, parser, intervals=args.intervals,
                                            auto_tabix=args.auto_tabix,
                            chunk=not args.no_chunk, chunk_factor=args.chunk_factor,
                                 tmp_dir=args.temp_dir, ncpus=args.ncpus,
                           skip_multiple_entries=args.skip_multiple_entries,
                            progress_report_interval=args.progress_report_interval)

    logging.info("Using {} to traverse the file.".format(walker.__class__.__name__))

    start = time.time()
    walker.run()
    end = time.time()
    delta = end - start
    logging.info("This run took {} seconds = {} minutes = {} hours.".format(delta,delta/60.,delta/3600.))



def main():
    argparser = get_argparser()
    args, sub_args = parse_args(argparser)
    try:
        parse(args,sub_args)
        return 0
    except KeyboardInterrupt:
        ### handle keyboard interrupt ###
        return 0
    except Exception, e:
        logging.exception(e)
        return 2


if __name__ == "__main__":
    sys.exit(main())








