#!/usr/bin/env python
"""
Kyle Hernandez
khernandez@bsd.uchicago.edu
Center for Research Informatics
University of Chicago
"""
import sys

def usage():
    """Prints usage to the screen"""
    print """
-------------------------------------------------------------------------------
Author: Kyle Hernandez <khernandez@bsd.uchicago.edu>

Description: Processes the biallelic, normalized RIL variants from freebayes 
and converts them to A/B calls based on the parental genotypes. Here, sample
Bd21 is called parent A and sample Bd3-1 is called parent B. 

Notes: 
- Only considers genotypes with a GQ >= 8
- Since Bd21 was also used to build the reference, missing parental data was
inferred where possible.

Usage:
	FreebayesToMarkers.py <freebayes.biallelic.vcf> <markers.raw.txt>
-------------------------------------------------------------------------------
"""

def getGenotype(d, g):
    curr = dict(zip(d, g.split(':')))
    if curr['GT'] == '.': return None
    elif int(curr['GQ']) < 8: return None
    else: return curr['GT']

def is_het(g):
    if g and g[0] != g[-1]: return True
    else: return False

def getMarkers(pa, pb, gts):
    curr = []
    if pa and pb:
        for g in gts:
            if g is None: curr.append('N')
            elif g == pa: curr.append('A')
            elif g == pb: curr.append('B')
            elif is_het(g): curr.append('H')
            else: curr.append('N')
    elif pa and pa == '0/0':
        for g in gts:
            if g is None: curr.append('N')
            elif g == pa: curr.append('A')
            elif is_het(g): curr.append('H')
            elif g == '1/1': curr.append('B')
            else: curr.append('N')
    elif pb and pb == '1/1':
        for g in gts:
            if g is None: curr.append('N')
            elif g == pb: curr.append('B')
            elif is_het(g): curr.append('H')
            elif g == '0/0': curr.append('A')
            else: curr.append('N')
    elif pa is None and pb is None:
        for g in gts:
            if g is None: curr.append('N')
            elif g == '0/0': curr.append('A')
            elif is_het(g): curr.append('H')
            elif g == '1/1': curr.append('B')
            else: curr.append('N')
    return curr


if __name__ == '__main__':
    if len(sys.argv) != 3:
        usage()
        sys.exit(1)

    # Command line args
    ifil = sys.argv[1]
    ofil = sys.argv[2]

    # Start processing file    
    with open(ofil, 'wb') as o:
        head = []
        samples = []

        # Parent Sample IDs
        parentA = 'BD_Bd21'
        parentB = 'BD_Bd3-1'

        # process input vcf
        for line in open(ifil, 'rU'):
            if line.startswith('##'): continue
            elif line.startswith('#CHROM'):
                head = line.rstrip().split('\t')
                samples = head[9:-2]
                o.write('\t'.join(['marker'] + samples) + '\n')
            else:
                dat = dict(zip(head, line.rstrip().split('\t')))
                fmt = dat['FORMAT'].split(':')
                pgt_A = getGenotype(fmt, dat[parentA])
                pgt_B = getGenotype(fmt, dat[parentB])
                gt_calls = [getGenotype(fmt, dat[i]) for i in samples]
                markers = []
                if is_het(pgt_A) or is_het(pgt_B): continue
                elif (pgt_A and pgt_B) and (pgt_A == pgt_B): continue
                else:
                    markers = getMarkers(pgt_A, pgt_B, gt_calls)
                    if markers:
                        o.write('\t'.join(['{#CHROM}:{POS}'.format(**dat)] + markers) + '\n')
