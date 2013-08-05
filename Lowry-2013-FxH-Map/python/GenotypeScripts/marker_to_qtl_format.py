#!/usr/bin/env python
# Kyle Hernandez
# marker_to_qtl_format.py - Filter the markers from the marker matrix created
#                           with vcf_to_marker.py and create matrix in various formats
#                           for R/qtl or MSTmap
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
import sys
import os
import random
import time

FMAT_QTL = {'mst_map': {'A': 'A', 'B': 'B', 'H': 'X', '-': 'U'},
            'join_map': {'A': 'a', 'B': 'b', 'H': 'h', '-': '-'}}

def main():
    '''
    ------------------------------------------------------------------------
    Author: Kyle Hernandez
    Email:  kmhernan@utexas.edu

    Take the marker matrix (output from vcf_to_marker.py) and filter sites
    based on % missing and segregation distortion. Print to various formats.
    Automatically filters to best hit when multiple hits are found on a RAD
    cutsite.
    ------------------------------------------------------------------------

    USAGE: marker_to_qtl_format.py matrix.tab out.file map_type max_missing seq_p reduce_scaff

    ARGUMENTS:
    	matrix.tab   - output matrix from vcf_to_marker.py
        out.file     - name of output file
        map_type     - Type of output file [mst_map | rqtl | join_map]
        max_missing  - Limit for % of data missing (Float)
        seg_p        - Limit for segregation distortion p-value (Float)
        reduce_scaff - (Y/N) If there are multiple hits/scaff, print out
                       the 'best' marker # NOT YET IMPLEMENTED
    '''
    header, site_dict = process_matrix()
    best_dict = reduce_RAD(site_dict)
    if map_type == 'mst_map':
        print_mst(header, best_dict)
    elif map_type == 'rqtl':
        print_rqtl(header, best_dict)
    elif map_type == 'join_map':
        print_join_map(header, best_dict)

def process_matrix():
    '''
       Build a dictionary of tuples containing 
       information about markers that pass filters.
       Return a list of the species columns and
       dictionary.
    '''
    header = []
    site_dict = {}

    with open(in_matrix, 'rU') as f:
        for line in f:

            if line.startswith('Loci'):
                header = line.rstrip().split('\t')[8::]

            else:
                cols      = line.rstrip().split('\t')
                sc,po     = [int(i) for i in cols[0].split(':')]
                fox       = cols[1]
                p_missing = float(cols[5])
                seg_p     = float(cols[7])
                value     = (fox, p_missing, seg_p, cols[8::])
 
                if p_missing <= max_missing:
                    if seg_p >= seg_limit:
                        if sc not in site_dict: site_dict[sc] = {}
                        site_dict[sc][po] = value 

    return header, site_dict

def reduce_RAD(site_dict):
    '''Reduce to only one snp/RAD tag and return dictionary of best sites'''
    best = {}

    for s in sorted(site_dict):
        curr_positions = sorted(site_dict[s].keys())

        if len(curr_positions) == 1:
            if s not in best: best[s] = {}
            best[s][curr_positions[0]] = site_dict[s][curr_positions[0]]

        else:
            if s not in best: best[s] = {}

            good = [i for i in curr_positions \
                   if len(set(range(i - 36, i + 37)) \
                   & set(curr_positions)) == 1]
            for i in good: best[s][i] = site_dict[s][i]

            multiple = list(set(curr_positions) - set(good))
            if multiple:
                if max(multiple) - min(multiple) <= 36:
                   L = sorted(multiple, key=lambda k: site_dict[s][k][1])
                   p_set = set([site_dict[s][i][1] for i in L[0:2]])
                   s_set = set([site_dict[s][i][2] for i in L[0:2]])

                   if len(p_set) == 1 and len(s_set) == 1:
                       # Make a random choice
                       choice = random.sample(L[0:2], 1)
                       best[s][choice[0]] = site_dict[s][choice[0]]

                   elif len(p_set) == 1 and len(s_set) > 1:
                       # Choose the one with the greatest seg dist pvalue
                       K = sorted(L, key=lambda k: site_dict[s][k][2])
                       choice = K.pop()
                       best[s][choice] = site_dict[s][choice]

                   else:
                       # Choose lowest p missing
                       choice = L.pop()
                       best[s][choice] = site_dict[s][choice]

                else:
                    curr_rad = []
                    for i in sorted(multiple):

                        if not curr_rad:
                            curr_rad.append(i)

                        elif curr_rad and i - curr_rad[-1] <= 36:
                            curr_rad.append(i)

                        else:
                           # Process 
                           L = sorted(curr_rad, key=lambda k: site_dict[s][k][1])
                           p_set = set([site_dict[s][j][1] for j in L[0:2]])
                           s_set = set([site_dict[s][j][2] for j in L[0:2]])
                      
                           if len(p_set) == 1 and len(s_set) == 1:
                    	       # Make a random choice
                       	       choice = random.sample(L[0:2], 1)
                               best[s][choice[0]] = site_dict[s][choice[0]]

                           elif len(p_set) == 1 and len(s_set) > 1:
                               # Choose the one with the greatest seg dist pvalue
                               K = sorted(L, key=lambda k: site_dict[s][k][2])
                               choice = K.pop()
                               best[s][choice] = site_dict[s][choice]

                           else:
                               # Choose lowest p missing
                               choice = L.pop()
                               best[s][choice] = site_dict[s][choice]
                           
                           # Reset for next rad tag 
                           curr_rad = []
                           curr_rad.append(i)

                    # Process last tag 
                    L = sorted(curr_rad, key=lambda k: site_dict[s][k][1])
                    p_set = set([site_dict[s][j][1] for j in L[0:2]])
                    s_set = set([site_dict[s][j][2] for j in L[0:2]])
                      
                    if len(p_set) == 1 and len(s_set) == 1:
                       # Make a random choice
                       choice = random.sample(L[0:2], 1)
                       best[s][choice[0]] = site_dict[s][choice[0]]

                    elif len(p_set) == 1 and len(s_set) > 1:
                       # Choose the one with the greatest seg dist pvalue
                       K = sorted(L, key=lambda k: site_dict[s][k][2])
                       choice = K.pop()
                       best[s][choice] = site_dict[s][choice]

                    else:
                       # Choose lowest p missing
                       choice = L.pop()
                       best[s][choice] = site_dict[s][choice]

    return best

def print_mst(header, dic):
    '''Prints out markers in mst_map format'''

    with open(out_matrix, 'wb') as o:
        print "Number of potential markers:", sum(len(x) for x in dic.values())
        hd = ['population_type RIL2','population_name HxF', 'distance_function kosambi',
              'cut_off_p_value 0.000001', 'no_map_dist 20', 'no_map_size 0', 'missing_threshold 0.8',
              'estimation_before_clustering no', 'detect_bad_data no', 'objective_function ML',
              'number_of_loci ' + str(sum(len(x) for x in dic.values())),
              'number_of_individual ' + str(len(header))]
        o.write('\n'.join(hd) + '\n' + 'Loci\t' + '\t'.join(header) + '\n')

        for s in sorted(dic):
            for p in sorted(dic[s]):
                mname = dic[s][p][0] + ':' + str(s) + ':' + str(p)
                o.write(mname + '\t' + '\t'.join([FMAT_QTL['mst_map'][i] for i in dic[s][p][3]]) + '\n')

def print_rqtl(header, dic):
    '''Prints out markers in rqtl csvr format [tab delim]'''

    with open(out_matrix, 'wb') as o:
        print "Number of potential markers:", sum(len(x) for x in dic.values())
        hd = ['Sample', ''] + header
        o.write('\t'.join(hd) + '\n')
            
        for s in sorted(dic):
            for p in sorted(dic[s]):
                mname = dic[s][p][0] + ':' + str(s) + ':' + str(p)
                o.write(mname+ '\t1\t' + '\t'.join(dic[s][p][3]) + '\n') 

def print_join_map(header, dic):
    '''Prints out markers in JoinMap loc-file format'''

    print "Number of potential markers:", sum(len(x) for x in dic.values())
    jmap_header = ['name = FxH', 'popt = F2', 
                   'nloc = ' + str(sum(len(x) for x in dic.values())),
                   'nind = ' + str(len(header))]

    with open(out_matrix, 'wb') as o:
        # output header
        o.write('\n'.join(jmap_header) + '\n\n')
        
        for s in sorted(dic):
            for p in sorted(dic[s]):
                mname = dic[s][p][0] + ':' + str(s) + ':' + str(p)
                o.write(mname + '\n')
                o.write('\t' + ''.join([FMAT_QTL['join_map'][i] for i in dic[s][p][3]]) + '\n')
         
        o.write('\nindividual names:\n\n')
        o.write('\n'.join(header) + '\n')

if __name__ == '__main__':
    fmat = ['mst_map', 'rqtl', 'join_map']
    if len(sys.argv) != 7 or sys.argv[3] not in fmat:
        print main.__doc__
        sys.exit()
    start        = time.time()
    in_matrix    = sys.argv[1]
    out_matrix   = sys.argv[2]
    map_type     = sys.argv[3]
    max_missing  = float(sys.argv[4])
    seg_limit    = float(sys.argv[5])
    reduce_scaff = sys.argv[6]
    main()
    print "Finished; Took:", time.time() - start, "seconds."
