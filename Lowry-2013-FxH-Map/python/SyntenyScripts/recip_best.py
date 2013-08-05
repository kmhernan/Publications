#!/usr/bin/env python
# Kyle Hernandez
#
# RecipBest.py - Parses the best hits.
# Although not required in any sense, share the love and pass on attribution
# when using or modifying this code.
#
# To the extent possible under law, the author(s) have dedicated all copyright 
# and related and neighboring rights to this software to the public domain 
# worldwide. This software is distributed without any warranty.
# 
# You should have received a copy of the CC0 Public Domain Dedication along with 
# this software. If not, see <http://creativecommons.org/publicdomain/zero/1.0/>
#
import time
import sys
from operator import itemgetter
import random
import collections

def main():
    """
    ---------------------------------------------------------------------------
    AUTHOR: Kyle Hernandez
    EMAIL: kmhernan@utexas.edu
   
    Parses out the reciprocal best hits from filtered Blast hits via ParseBlast.py.
    ---------------------------------------------------------------------------

    USAGE: RecipBest.py fh_fox_best.tab fox_fh_best.tab out.tab

    ARGUMENTS:
         fh_fox_best.tab - tab delimited best hits file of fh blasted onto Foxtail 
         fox_fh_best.tab - tab delimited best hits file of fox blasted onto FH
         out.tab   - tab delimited best hits file
    """
    FH_dict = build_FH()
    process_Fox(FH_dict)

def build_FH():
    """
    Creates the dictionary of FH best hits on foxtail 
    """
    dic = {} 
    with open(fh_fox_fil, 'rU') as f:
        for line in f:
            cols   = line.rstrip().split('\t')
            query  = cols[0]
            target = cols[1]

            if query not in dic: dic[query] = {}
            dic[query][target] = ''
    return dic

def process_Fox(dic):
    """
    Finds the best recip hits and writes to file 
    """
    num_match = 0
    b_dic = {}
    o = open(out_fil, 'wb')

    with open(fox_fh_fil, 'rU') as f:
        for line in f:
            cols   = line.rstrip().split('\t')
            query  = cols[0]
            target = cols[1]
            evalue = float(cols[10])

            ch = query.split('_')[0]
            if target in dic:
                if ch in dic[target]:
                    num_match += 1
                    tvals = target.split('_')
                    o.write(ch + '\t' + '\t'.join(tvals) + '\n') 

    o.close()    
    print "Number of best hits:", num_match
 
if __name__ == '__main__':
    start = time.time()
    if len(sys.argv) != 4:
        print main.__doc__
        sys.exit()
    fh_fox_fil = sys.argv[1]
    fox_fh_fil = sys.argv[2]
    out_fil    = sys.argv[3]
    main()
    print "Finished; Took:", time.time() - start, "seconds."
