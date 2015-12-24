#!/usr/bin/env python
"""
Kyle Hernandez
khernandez@bsd.uchicago.edu
Center for Research Informatics
University of Chicago
"""
import sys
from scipy import stats 

def usage():
    """Prints usage to the screen"""
    print """
-------------------------------------------------------------------------------
Author: Kyle Hernandez <khernandez@bsd.uchicago.edu>

Description: Gets the top 90th percentile high-quality variants from a VCF file
and outputs them to a new VCF file. This method is used for the base quality
recalibration methods 

Notes: 
- Only processes loci with 'PASS' in the FILTER column
- Only processes loci where at least one sample has a SNP GT > 90
Usage:
	GetHighQualSNPs.py <input.vcf> <output.vcf> 
-------------------------------------------------------------------------------
"""

def extract_quals():
    """Extracts the quality scores and returns the 90th percentile cutoff"""
    AQ = []
    for line in open(fil, 'rU'):
        if line.startswith('#'): continue
        else:
            cols = line.rstrip().split('\t')
            flt  = cols[6]
            if flt == 'PASS':
                AQ.append(float(cols[5])) 
    return stats.scoreatpercentile(AQ, 90)

def filter_vcf(percentile):
    """Applies the filtering"""
    with open(ofil, 'wb') as o:
        head = []
        samples = []
        for line in open(fil, 'rU'):
            if line.startswith('##'): o.write(line) 
            elif line.startswith('#'):
                o.write(line)
                head = line.rstrip().split('\t')
                samples = head[9:]
            else:
                cols = line.rstrip().split('\t')
                flt  = cols[6]
                qual = float(cols[5])
                if flt == 'PASS' and qual > percentile:
                    fmt = cols[8].split(':')
                    idx = fmt.index('GQ')
                    called = [i.split(':') for i in cols[9:] if i.split(':')[0] != './.']
                    good   = [i for i in called if i[0] != '0/0' and int(i[idx]) > 90] 
                    if good: 
                        o.write(line)

if __name__ == '__main__':
    if len(sys.argv) != 3:
        usage()
        sys.exit(1)

    # Command line args
    fil  = sys.argv[1] 
    ofil = sys.argv[2] 

    # Get 90th percentile cutoff
    percentile = extract_quals()

    # Apply filters
    filter_vcf(percentile)
