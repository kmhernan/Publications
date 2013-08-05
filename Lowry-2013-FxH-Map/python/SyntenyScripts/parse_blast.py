#!/usr/bin/env python
# Kyle Hernandez
#
# ParseBlast.py - Parses the best hits.
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
   
    Parses out the best hits from Blastn
    ---------------------------------------------------------------------------

    USAGE: ParseBlast.py blast.tab out.tab

    ARGUMENTS:
         blast.tab - tab delimited Blast output file
         out.tab   - tab delimited best hits file
    """
    blast_dict = blastdict_builder()
    process_dict(blast_dict)

def blastdict_builder():
    """
    Creates the dictionary of blast hits
    """
    dic = {} 
    with open(blast_fil, 'rU') as f:
        for line in f:
            cols   = line.rstrip().split('\t')
            query  = cols[0]
            target = cols[1]
            evalue = float(cols[10])

            if evalue <= float(10e-10):
               data = [int(i) if n == 6 or n == 7 else float(i) if n == 8 or n == 9 else i for n,i in enumerate(cols[2::])]
               if query not in dic: dic[query] = {}
               if target not in dic[query]: dic[query][target] = []
               dic[query][target].append(data)
    return dic

def process_dict(dic):
    """
    Finds the best hit based on the filtered hits in the dictionary.
    """
    num_match = 0
    with open(out_fil, 'wb') as o:
        for q in sorted(dic):
            chr_hits = sorted(dic[q].keys())

            if len(chr_hits) == 1:
                curr_data = dic[q][chr_hits[0]]

                if len(curr_data) == 1:
                    num_match += 1
                    o.write(q + '\t' + chr_hits[0] + '\t' + \
                           '\t'.join([str(i) for i in curr_data[0]]) + '\n')
    
                else:
                    L = sorted(curr_data, key=itemgetter(8))
                    
                    if len(set([i[8] for i in L[0:2]])) == 2:
                        num_match += 1
                        o.write(q + '\t' + chr_hits[0] + '\t' + \
                           '\t'.join([str(i) for i in L[0]]) + '\n')
                    elif L[0][9] > L[1][9]:
                        num_match += 1
                        o.write(q + '\t' + chr_hits[0] + '\t' + \
                           '\t'.join([str(i) for i in L[0]]) + '\n')
                    elif L[0][9] < L[1][9]:
                        num_match += 1
                        o.write(q + '\t' + chr_hits[0] + '\t' + \
                           '\t'.join([str(i) for i in L[0]]) + '\n')
                    elif len(set([i[8] for i in L])) == 1: 
                        choice = random.sample(L,1).pop()
                        num_match += 1
                        o.write(q + '\t' + chr_hits[0] + '\t' + \
                           '\t'.join([str(i) for i in choice]) + '\n')
                    else: continue
            # Deal with sequences that blast to multiple scaffolds
            # finding the scaffold with the minimum sum e-value
            else:
                curr_sums = {} 
                for c in chr_hits:
                    if c not in curr_sums: 
                        curr_sums[c] = 0
                    curr_data = dic[q][c]
                    curr_sums[c] += sum([i[8] for i in curr_data])
                
                L = sorted(curr_sums.iteritems(), key=itemgetter(1))
                
                if L[0][1] < L[1][1]:
                    best_ch = L[0][0]
                    
                    best_scaff = dic[q][best_ch]
                    
                    if len(best_scaff) == 1:
                        
                        num_match += 1
                        o.write(q + '\t' + best_ch + '\t' + \
                               '\t'.join([str(i) for i in best_scaff[0]]) + '\n')
                    
                    else:
                        K = sorted(best_scaff, key=itemgetter(8))
                    
                        if len(set([i[8] for i in K[0:2]])) == 2:
                            num_match += 1
                            o.write(q + '\t' + L[0][0] + '\t' + \
                               '\t'.join([str(i) for i in K[0]]) + '\n')

                        elif K[0][9] > K[1][9]:
                            num_match += 1
                            o.write(q + '\t' + L[0][0] + '\t' + \
                                '\t'.join([str(i) for i in K[0]]) + '\n')
                    
                        elif K[0][9] < K[1][9]:
                            num_match += 1
                            o.write(q + '\t' + L[0][0] + '\t' + \
                                '\t'.join([str(i) for i in K[0]]) + '\n')
                   
                        elif len(set([i[8] for i in K])) == 1: 
                            choice = random.sample(K,1).pop()
                            num_match += 1
                            o.write(q + '\t' + L[0][0] + '\t' + \
                                '\t'.join([str(i) for i in choice]) + '\n')
            
                        else: continue
    
    print "Number of best hits:", num_match
 
if __name__ == '__main__':
    start = time.time()
    if len(sys.argv) != 3:
        print main.__doc__
        sys.exit()
    blast_fil = sys.argv[1]
    out_fil   = sys.argv[2]
    main()
    print "Finished; Took:", time.time() - start, "seconds."
