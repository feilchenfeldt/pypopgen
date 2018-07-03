#!/usr/bin/env python
"""
Reads a line from 'samtools mpileup -s'
and prints chrom, pos, rint coverage, 
root mean square
mapping quality and number of
mapping quality zero reads for
a mapping quality string in
phred33 encoding.
"""
import sys, math

l = sys.argv[1].strip().split()
chrom = l[0]
pos = l[1]

if l[0] == '*':
    nmq0 = '*'
    rmsmq = '*'
    cov = 0
else:
    mqs = [ord(s)-33 for s in l[6]]
    nmq0 = mqs.count(0)
    mqs1 = [m for m in mqs if m!=0]
    rmsmq = math.sqrt(sum([m**2 for m in mqs1])*1./len(mqs1))
    cov = len(mqs)

print '\t'.join([chrom, pos, str(cov), str(rmsmq), str(nmq0)])


