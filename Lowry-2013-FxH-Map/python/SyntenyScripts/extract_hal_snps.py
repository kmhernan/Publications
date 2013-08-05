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
import re

def main():
    """
    ---------------------------------------------------------------------------
    AUTHOR: Kyle Hernandez
    EMAIL: kmhernan@utexas.edu
    
    Locates the SNPs falling in best hits exons of FTM. 
    ---------------------------------------------------------------------------

    USAGE: extract_hal_snps.py inbest.tab exon.fa hxf.vcf out.tab 

    ARGUMENTS:
    	inbest.tab  - Tab delimited output from Blast Parser of Fox exons on Hallii
      	exon.fa     - Fasta file of exons with exon ID in > line
        hxf.vcf     - HxF VCF file.
        out.tab     - The output file
    """
    exon_data = get_exon_data()
    cut_dict = get_cut_sites(exon_data)
    find_snps(cut_dict)

def get_exon_data():
    """
    Loads exon data to dic
    """
    dic = {}
    with open(exonfa, 'rU') as f:
        for line in f:
            if line.startswith('>'):
                cols = line.rstrip().split(' | ')
                exon = cols[0][1::]
                dat = re.split( r'pacid=', cols[1])[1]
                dic[exon] = dat
    return dic

def get_cut_sites(exon):
    """
    Reads the blast output and returns a dict of start/stop positions 1-index scaled
    """
    dic = {}

    with open(inbest, 'rU') as f:
        for line in f:
            cols = line.rstrip().split('\t')
            query = cols[0]
            ch   = cols[1]
            p1   = min([int(cols[8]), int(cols[9])])
            p2   = max([int(cols[8]), int(cols[9])])
            
            val  = (p1, p2, query, exon[query])
            if ch not in dic: dic[ch] = {} 
            dic[ch][val] = [cols[10],cols[11]]

    return dic

def find_snps(cut_dict):
    """
    Finds SNPs within the VCF file that occur in the best hit blast regions.
    """
    dic = {}
    o = open(outfil, 'wb')
    o.write('FHCHRM\tFHPOS\tFM_PACID\tFMCHRM\tFMORI\tFMSTART\tFMEND\tEVALUE\tSCORE\tKEY\n')
    with open(invcf, 'rU') as f:
        for line in f:
            if line.startswith('#'): continue
            cols = line.rstrip().split('\t')
            if cols[4] != '.':
                ch   = cols[0]
                po   = int(cols[1])
           
                try: 
                    curr_cut = cut_dict[ch]
                    L = [i for i in curr_cut.keys() if i[0] <= po <= i[1]]
                    if L:
                        for j in L:
                            query = j[2].split('_')
                            fmchrm = query[1]
                            fmst   = query[2]
                            fmen   = query[3]
                            fmori  = query[4]
                            row = [ch, str(po), j[3], fmchrm, fmori, str(j[0]), str(j[1]), 
                                   cut_dict[ch][j][0], cut_dict[ch][j][1], j[2]]
                            o.write('\t'.join(row) + '\n')
                except KeyError: pass 
    o.close()

if __name__ == '__main__':
    start = time.time()
    if len(sys.argv) != 5:
        print main.__doc__
        sys.exit()
    inbest  = sys.argv[1]
    exonfa  = sys.argv[2]
    invcf   = sys.argv[3]
    outfil  = sys.argv[4]
    main()
    print "Finished; Took:", time.time() - start, "seconds."
