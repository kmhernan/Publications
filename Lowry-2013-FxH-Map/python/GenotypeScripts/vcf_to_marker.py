#!/usr/bin/env python
# Kyle Hernandez, 2013
# kmhernan@utexas.edu 
# vcf_to_marker.py - Converts the multi-sample VCF file of the FH mapping population
#                    into markers for QTL mapping.
#
# WARNING!! Not generalizable, specific to our FH mapping sequences
#
# This is free and unencumbered software released into the public domain.
# 
# Anyone is free to copy, modify, publish, use, compile, sell, or
# distribute this software, either in source code form or as a compiled
# binary, for any purpose, commercial or non-commercial, and by any
# means.
# 
# In jurisdictions that recognize copyright laws, the author or authors
# of this software dedicate any and all copyright interest in the
# software to the public domain. We make this dedication for the benefit
# of the public at large and to the detriment of our heirs and
# successors. We intend this dedication to be an overt act of
# relinquishment in perpetuity of all present and future rights to this
# software under copyright law.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
# IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR
# OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
# ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
# OTHER DEALINGS IN THE SOFTWARE.
# 
# For more information, please refer to <http://unlicense.org/>

import sys
from scipy import stats
import numpy as np
import time

# Global Parent ID Variables
par_A = 'Hal2'
par_B = 'FIL2'
F1    = 'FHF1'

# Global Collections
CHRM_DICT = {}
ct_dict = {
  'Total SNPs': 0,
  'Not <aaxbb> Marker': 0,
  'Potential Marker': 0,
  'Low GQ': 0,
  'Not SNP': 0,
  'Not biallelic': 0
}

def main():
    '''
    --------------------------------------------------------
    Author:  Kyle Hernandez
    Contact: kmhernan@utexas.edu

    Take the multi-sample mapping population VCF file,
    find biallelic SNPs, remove low GQs, call genotypes
    relative to Hal (A) and Fil (B), estimate various 
    site statistics, add in Foxtail Millet BLAST information
    and output to tab delimited flat file.
    ---------------------------------------------------------

    USAGE: vcf_to_marker.py <file.vcf> <BLAST.tab> <out.tab> <GQ_limit>

    ARGUMENTS:
    	file.vcf  - VCF file from FH reads
        BLAST.tab - Formatted recip best-hits output from RecipBest
        out.tab   - output file 
        GQ_limit  - Minimally acceptable GQ (Int)
    '''
    load_foxtail()
    process_vcf()
    print_summary()

def load_foxtail():
    '''Loads BLAST information'''
    with open(blast_file, 'rU') as f:
        for line in f:
            cols = line.rstrip().split('\t')
            ch = cols[0].replace('chr','')
            scf = cols[1]
            key = (int(cols[2]), int(cols[3]))

            if scf not in CHRM_DICT: CHRM_DICT[scf] = {}
            CHRM_DICT[scf][key] = ch

def process_vcf():
    '''Process the VCF file'''
    o = open(out_file, 'wb')

    with open(vcf_file, 'rU') as f:
        for line in f:
            # Skip meta-information
            if line.startswith('##'): continue
            # Grab header line
            elif line.startswith('#CHROM'):
                samples_list, kids_list = get_samples(line)
                o.write('Loci\tFoxtail_Chr\tMarker_Class\tAvg_DP\tAvg_GQ\tP_Missing\tA:H:B\tSeg_Pval\t' + \
                        '\t'.join(kids_list) + '\n')

            else:
                # initialize containers for each SNP
                base = {}
                cols = line.rstrip().split('\t')

                # Apply all filters and count
                if cols[4] != '.':
                    if ',' not in cols[4]:
                        # Process biallelic SNP
                        ct_dict['Total SNPs'] += 1
                        base = process_bases(samples_list, base, cols)
                        parent_pattern = base[par_A]['gt'] +' '+ base[par_B]['gt']
                 
                        # Get marker class and marker calls
                        f2_genotypes = {}
                        marker_class = ''
                        marker_class,\
                        f2_genotypes = get_marker(marker_class, parent_pattern, f2_genotypes)
                 
                        # Process good markers                   
                        if marker_class:
                            loci = cols[0] + ':' + cols[1]
                            print_row_with_stats(marker_class, f2_genotypes, base, kids_list, loci, o)

                        else:
                            ct_dict['Not <aaxbb> Marker'] += 1
                    else:
                        ct_dict['Not biallelic'] += 1
                else:
                    ct_dict['Not SNP'] += 1
    o.close()

def get_samples(line):
    '''Get header'''
    samples_list = line.rstrip().split('\t')[9::]
    kids_list    = [i for i in samples_list if i != par_A \
                    and i != par_B and i != F1]
    return samples_list, kids_list

def process_bases(samples_list, base, cols):
    '''Process current sample and return dict of F2 GTs, depth, and quality'''
    for i,s in enumerate(samples_list):
        # Deal with samples that are Ns
        if cols[9::][i] == './.':
            if s not in base: base[s] = {}
            base[s]['gt'] = '-'
            base[s]['dp'] = '-'
            base[s]['gq'] = '-'

        else:
            sample_fields = cols[9::][i].split(':')
             
            # Does it pass GQ minimum?
            if int(sample_fields[3]) >= GQ_limit: 
                if s not in base: base[s] = {}
                base[s]['gt'] = sample_fields[0]
                base[s]['dp'] = int(sample_fields[2])
                base[s]['gq'] = int(sample_fields[3])

            else:
                ct_dict['Low GQ']        += 1
                if s not in base: base[s] = {}
                base[s]['gt'] = '-'
                base[s]['dp'] = '-'
                base[s]['gq'] = '-'
    return base

def get_marker(marker_class, parent_pattern, f2_genotypes):
    ''' 
        I currently only deal with the simple AAxBB markers. We
        can potentially use other markers like:
        0/0 0/1, 1/1 0/1, 0/1 0/0, 0/1 1/1, 0/1 0/1 and possibly
        sites with more than 2 alleles if we use the F1 data
     '''
    if parent_pattern == '0/0 1/1':
        marker_class        = '<aaxbb>'
        f2_genotypes['0/0'] = 'A'
        f2_genotypes['0/1'] = 'H'
        f2_genotypes['1/0'] = 'H'
        f2_genotypes['1/1'] = 'B'
        f2_genotypes['-']   = '-'
 
    elif parent_pattern == '1/1 0/0':
        marker_class        = '<aaxbb>'
        f2_genotypes['0/0'] = 'B'
        f2_genotypes['0/1'] = 'H'
        f2_genotypes['1/0'] = 'H'
        f2_genotypes['1/1'] = 'A'
        f2_genotypes['-']   = '-'

    return marker_class, f2_genotypes

def print_row_with_stats(marker_class, f2_genotypes, base, kids_list, loci, o):
    '''Print out marker stats and genotypes'''
    curr_GTs = [f2_genotypes[base[i]['gt']] for i in kids_list]

    if len(set(curr_GTs)) > 1:

        ct_dict['Potential Marker'] += 1

        # Some site statistics
        avg_depth = np.mean(np.array([base[i]['dp']\
                                      for i in kids_list\
                                      if base[i]['dp'] != '-']))
        avg_gq    = np.mean(np.array([base[i]['gq']\
                                      for i in kids_list\
                                      if base[i]['gq'] != '-']))
        perc_missing   = float(curr_GTs.count('-'))/len(kids_list)
        observed, pval = get_seg_dist(curr_GTs)

        # Get Foxtail information
        try:
            scf = loci.split(':')[0]
            st  = int(loci.split(':')[1])
            curr_pos = CHRM_DICT[scf]
            
            fox_c = 'NA'
            for i in sorted(curr_pos.keys()):
                if st in range(i[0], i[1]):
                    fox_c = curr_pos[i] 
        except:
            fox_c = 'NA'
        # Write to flat file
        o.write('\t'.join([loci, fox_c, marker_class, 
                    '{:0.2F}'.format(avg_depth), 
                    '{:0.2F}'.format(avg_gq),
                    '{:0.2F}'.format(perc_missing),
                    observed, '{:0.2E}'.format(pval)] +\
                    curr_GTs) + '\n')

def get_seg_dist(GTs):
    '''Checks for segregation distortion. Returns ratio string and pvalue'''
    observed = np.hstack([GTs.count('A'), GTs.count('H'), GTs.count('B')])
    n_typed  = len(GTs) - GTs.count('-')
    expected = np.hstack([0.25, 0.50, 0.25])
    ch2, pval = stats.chisquare(observed, expected*n_typed)
    obs = ':'.join([str(round(i/float(n_typed), 2)) for i in observed])
    return obs, pval

def print_summary():
    print 'Monomorphic:', ct_dict['Not SNP']
    print 'Total SNPs:', ct_dict['Total SNPs']
    print 'Number of SNPs > 2 alleles:', ct_dict['Not biallelic']
    print 'Number of low GQ at threshold', GQ_limit, ':', ct_dict['Low GQ']
    print 'Not aaxbb Marker:', ct_dict['Not <aaxbb> Marker']    
    print 'Potential Markers:', ct_dict['Potential Marker']

if __name__=='__main__':
    if len(sys.argv) != 5:
        print main.__doc__
        sys.exit()
    start      = time.time()
    vcf_file   = sys.argv[1]
    blast_file = sys.argv[2]
    out_file   = sys.argv[3]
    GQ_limit   = int(sys.argv[4])
    main()
    print "Finished; Took:", time.time() - start, "seconds."
