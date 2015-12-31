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

Description: Filters out the raw marker calls based on the frequency of missing
data and the allele frequencies. 

Notes:
- Missing frequency cutoff of 0.30 
- A freq/B freq cutoff 0.60
- H freq cutoff 0.30 

Usage:
        FilterMarkers.py <markers.raw.txt> <markers.filtered.txt>
-------------------------------------------------------------------------------
"""

if __name__ == '__main__':
    if len(sys.argv) != 3:
        usage()
        sys.exit(1)

    # Get command line args
    ifil = sys.argv[1]
    ofil = sys.argv[2]

    # Process file 
    head = []
    with open(ofil, 'wb') as o:
        for line in open(ifil, 'rU'):
            if not head:
                head = line.rstrip().split('\t')
                o.write(line)
            else:
                cols  = line.rstrip().split('\t')
                calls = cols[1:]
                # Get missing freq
                fmiss = calls.count('N') / float(len(calls))
                if fmiss < 0.30:
                    # Get allele freqs
                    nomiss = [i for i in calls if i != 'N']
                    afreq  = nomiss.count('A') / float(len(nomiss))
                    bfreq  = nomiss.count('B') / float(len(nomiss))
                    hfreq  = nomiss.count('H') / float(len(nomiss))
                    if afreq < 0.60 and bfreq < 0.60 and hfreq < 0.30:
                        o.write(line)
