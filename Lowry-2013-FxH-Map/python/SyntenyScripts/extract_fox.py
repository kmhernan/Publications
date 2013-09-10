#!/usr/bin/env python
# Kyle Hernandez
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
import time

def main():
    """
    ---------------------------------------------------------------------------
    AUTHOR: Kyle Hernandez
    EMAIL: kmhernan@utexas.edu
 
    Takes the best hits of FH on Foxtail and extracts those locations from the
    Foxtail genome. Creates a fasta file of these sites for blasting.
    ---------------------------------------------------------------------------

    USAGE: ExtractFox.py inbest.tab infasta.fa outfasta.fa

    ARGUMENTS:
    	inbest.tab  - Tab delimited output from Blast Parser
        infasta.fa  - The Foxtail reference
        outfasta.fa - The output file for the cut Foxtail reference.
    """ 
    cut_dict = get_cut_sites()
    process_reference(cut_dict)

def get_cut_sites():
    """
    Reads the blast output and returns a dict of start/stop positions 0-index scaled
    """
    dic = {}

    with open(inbest, 'rU') as f:
        for line in f:
            cols = line.rstrip().split('\t')
            ch   = cols[1]
            p1   = min([int(cols[8]), int(cols[9])])
            p2   = max([int(cols[8]), int(cols[9])])
            if p1 - 1 < 0: start = 0
            else: start = p1 - 1
            
            val  = (start, p2)
            if ch not in dic: dic[ch] = []
            dic[ch].append(val)

    return dic

def process_reference(dic):
    """
    Reads in the Foxtail fasta reference file and cuts at given positions.
    """
    print 'Processing reference...'

    curr_scaff = ''
    curr_seq   = []
    o = open(outfasta, 'wb')

    with open(infasta, 'rU') as f:
        for line in f:
            if line.startswith('>') and not curr_scaff:
                curr_scaff = line.rstrip().split('>')[1]

            elif line.startswith('>') and curr_scaff:
                seq = ''.join(curr_seq)
                cuts = process_scaff(curr_scaff, dic[curr_scaff], seq)
                [print_cuts(j, curr_scaff, o) for j in cuts]

                curr_scaff = line.rstrip().split('>')[1]
                curr_seq   = []

            else:
                curr_seq.append(line.rstrip())
        seq = ''.join(curr_seq)
        cuts = process_scaff(curr_scaff, dic[curr_scaff], seq)
        [print_cuts(j, curr_scaff, o) for j in cuts]
    o.close()

def process_scaff(scf, pos, seq):
    """
    Returns cut sites.
    """
    for i in pos:
        if not seq[i[0]:i[1]]:
            print scf, i, len(seq)
    return [(i, seq[i[0]:i[1]]) for i in pos if seq[i[0]:i[1]]]

def print_cuts(c, s, o):
    """
    Write cut reference to file
    """
    o.write('>' + str(s) + '_' + str(c[0][0]) + '_' + str(c[0][1]) + '\n')
    o.write(c[1] + '\n')

if __name__ == '__main__':
    start = time.time()
    if len(sys.argv) != 4:
        print main.__doc__
        sys.exit()
    inbest   = sys.argv[1]
    infasta  = sys.argv[2]
    outfasta = sys.argv[3]
    main()
    print "Finished; Took:", time.time() - start, "seconds."
