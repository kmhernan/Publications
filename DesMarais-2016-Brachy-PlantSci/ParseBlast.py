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

Description: Gets the top hits from the Cui markers and converts them to a bed
file. 

Notes: This is a very simplistic script to only pull out the queries with a 
single hit.

Usage:
        ParseBlast.py <blast_output> <cui_markers.v3.0.bed>
-------------------------------------------------------------------------------
"""

def process_hit(q, d, o):
    '''
    Processes the current hit and writes to bed file
    '''
    if len(d) == 1:
        hit = d[0]
        row = [hit[1], str(min(int(hit[8]), int(hit[9]))), str(max(int(hit[8]), int(hit[9]))), q]
        o.write('\t'.join(row) + '\n')

if __name__ == '__main__':
    if len(sys.argv) != 3:
        usage()
        exit(1)

    # Parse command line
    infil = sys.argv[1]
    ofil  = sys.argv[2]

    # Process
    with open(ofil, 'wb') as o:

        # Set up record containers
        query = ''
        data = []

        # Process hits
        for line in open(infil, 'rU'):

            # First query 
            if line.startswith('# Query:') and not query:
                query = line.rstrip().split(' ')[-1]

            # All other queries 
            elif line.startswith('# Query:') and query:
                # Process hits from the query
                process_hit(query, data, o)

                # Reset to current query
                query = line.rstrip().split(' ')[-1]
                data  = []

            # Comment lines 
            elif line.startswith('#'): continue

            # Add the queries hit 
            else:
                data.append(line.rstrip().split('\t'))

        # Process last query
        process_hit(query, data, o)
