#!/usr/bin/env python
import sys
import time
from operator import itemgetter
from collections import OrderedDict
import reference

def main():
    """
    ---------------------------------------------------------------------------
    AUTHOR: Kyle Hernandez
    EMAIL: kmhernan@utexas.edu

    Call synteny by flanking blasts.
    ---------------------------------------------------------------------------
    USAGE: python flanking_exons.py best.tab link.tab ref.fa out.tab
    """
    best_dict = load_best()
    pt_dict = process_linkage(best_dict)
    print_data(pt_dict)

def load_best():
    """
    loads the best hits into a dict.
    """
    dic_fh = {}
    for line in open(inbest, 'rU'):
        if line.startswith('FHSCAF'): continue
        cols = line.rstrip().split('\t')
        fhch = int(cols[0])
        fhpo = int(cols[2])
        fmch = int(cols[1])
        fmst = int(cols[4])
        fmen = int(cols[5])
        evalue = float(cols[6])
        ftmkey = cols[9] 
        value  = (fmch, fmst, fmen, evalue, ftmkey)

        if fhch not in dic_fh: dic_fh[fhch] = {}
        if fhpo not in dic_fh[fhch]: dic_fh[fhch][fhpo] = [] 
        dic_fh[fhch][fhpo].append(value)
    return dic_fh

def process_linkage(best):
    """
    Loads linkage map data to dict
    """
    dic = OrderedDict() 

    for line in open(linkin, 'rU'):
        if line.startswith('Marker'): continue
        cols = line.rstrip().split('\t')
        ch, po  = map(int,cols[0].split(':'))
        lg = int(cols[1])
        cm = float(cols[2])
        
        mkr = cols[0]
        best_below = best_above = ''
        flank_below = flank_above = ftm_ch = ftm_key = 'NA' 
        try:
            curr_fm = sorted(best[ch][po])
            bst = ''
            if len(curr_fm) > 1:
                C = sorted(curr_fm, key=itemgetter(3))
                if C[0][3] < C[1][3]:
                    bst = C[0]
                elif C[0] == C[1]:
                    bst = C[0]
            else:
                bst = curr_fm[0]

            if bst:
                ftm_ch = bst[0]
                flank_below = bst[1]
                flank_above = bst[2]
                ftm_key = bst[4]
        except KeyError:
            try:
                curr_fm = sorted(best[ch].keys())
                try: below = min([i for i in curr_fm if i < po], key=lambda x:abs(x-po))
                except: below = None
                try: above = min([i for i in curr_fm if i > po], key=lambda x:abs(x-po)) 
                except: above = None
                
                if below:
                    below_dat = sorted(best[ch][below])

                    if len(below_dat) > 1:
                        B = sorted(below_dat, key=itemgetter(3))
                        if B[0][3] < B[1][3]:
                            best_below = B[0]
                        elif B[0] == B[1]:
                            best_below = B[0]
                    else:
                        best_below = below_dat[0]

                if above:
                    above_dat = sorted(best[ch][above])

                    if len(above_dat) > 1:
                        A = sorted(above_dat, key=itemgetter(3))
                        if A[0][3] < A[1][3]:
                            best_above = A[0]
                        elif A[0] == A[1]:
                            best_above = A[0]
                    else:
                        best_above = above_dat[0]

                if best_below and best_above:
                    if best_below[0] == best_above[0]:
                        ftm_ch = best_below[0]
                        flank_below = min([best_below[1], best_above[1]])
                        flank_above = max([best_below[2], best_above[2]])
                        ftm_key = best_below[4]
                    else:
                        AB = sorted([best_below, best_above], key=itemgetter(3))
                        if AB[0][3] < AB[1][3]:
                            ftm_ch = AB[0][0]
                            flank_below = min([AB[0][1], AB[0][1]])
                            flank_above = max([AB[0][2], AB[0][2]])
                            ftm_key = AB[0][4]
                elif best_below and not best_above:
                    ftm_ch = best_below[0]
                    flank_below = best_below[1]
                    flank_above = best_below[2]
                    ftm_key = best_below[4]
                elif best_above and not best_below:
                    ftm_ch = best_above[0]
                    flank_below = best_above[1]
                    flank_above = best_above[2]
                    ftm_key = best_above[4]
            except KeyError:
                pass
        if lg not in dic: dic[lg] = OrderedDict()
        if cm not in dic[lg]: dic[lg][cm] = OrderedDict()
        dic[lg][cm][mkr] = (lg, cm, mkr, ftm_ch, flank_below, flank_above, ftm_key)
    return dic
         
def process_linkage_old(best):
    """
    Loads linkage map data to dict
    """
    dic = OrderedDict() 

    for line in open(linkin, 'rU'):
        if line.startswith('Linkage'): continue
        cols = line.rstrip().split('\t')
        lg = int(cols[0])
        cm = float(cols[1])
        ch, po  = map(int,cols[2].split(':'))
        mkr = cols[2]
        best_below = best_above = ''
        flank_below = flank_above = ftm_ch = ftm_key = 'NA' 
        try:
            curr_fm = sorted(best[ch][po])
            bst = ''
            if len(curr_fm) > 1:
                C = sorted(curr_fm, key=itemgetter(3))
                if C[0][3] < C[1][3]:
                    bst = C[0]
                elif C[0] == C[1]:
                    bst = C[0]
            else:
                bst = curr_fm[0]

            if bst:
                ftm_ch = bst[0]
                flank_below = bst[1]
                flank_above = bst[2]
                ftm_key = bst[4]
        except KeyError:
            try:
                curr_fm = sorted(best[ch].keys())
                try: below = min([i for i in curr_fm if i < po], key=lambda x:abs(x-po))
                except: below = None
                try: above = min([i for i in curr_fm if i > po], key=lambda x:abs(x-po)) 
                except: above = None
                
                if below:
                    below_dat = sorted(best[ch][below])

                    if len(below_dat) > 1:
                        B = sorted(below_dat, key=itemgetter(3))
                        if B[0][3] < B[1][3]:
                            best_below = B[0]
                        elif B[0] == B[1]:
                            best_below = B[0]
                    else:
                        best_below = below_dat[0]

                if above:
                    above_dat = sorted(best[ch][above])

                    if len(above_dat) > 1:
                        A = sorted(above_dat, key=itemgetter(3))
                        if A[0][3] < A[1][3]:
                            best_above = A[0]
                        elif A[0] == A[1]:
                            best_above = A[0]
                    else:
                        best_above = above_dat[0]

                if best_below and best_above:
                    if best_below[0] == best_above[0]:
                        ftm_ch = best_below[0]
                        flank_below = min([best_below[1], best_above[1]])
                        flank_above = max([best_below[2], best_above[2]])
                        ftm_key = best_below[4]
                    else:
                        AB = sorted([best_below, best_above], key=itemgetter(3))
                        if AB[0][3] < AB[1][3]:
                            ftm_ch = AB[0][0]
                            flank_below = min([AB[0][1], AB[0][1]])
                            flank_above = max([AB[0][2], AB[0][2]])
                            ftm_key = AB[0][4]
                elif best_below and not best_above:
                    ftm_ch = best_below[0]
                    flank_below = best_below[1]
                    flank_above = best_below[2]
                    ftm_key = best_below[4]
                elif best_above and not best_below:
                    ftm_ch = best_above[0]
                    flank_below = best_above[1]
                    flank_above = best_above[2]
                    ftm_key = best_above[4]
            except KeyError:
                pass
        if lg not in dic: dic[lg] = OrderedDict()
        if cm not in dic[lg]: dic[lg][cm] = OrderedDict()
        dic[lg][cm][mkr] = (lg, cm, mkr, ftm_ch, flank_below, flank_above, ftm_key)
    return dic

def print_data(data):
    """
    Format for R plot and print.
    """
    ftm_ln = load_reflen() 
    o = open(outfil, 'wb')
    o.write('LG\tCm\tMarker\tFTM_Ch\tFTM_P1\tFTM_P2\trcm\trp1\trp2\n')
    for lg in data:
        mxcm = max(data[lg].keys())
        for cm in data[lg]:
            rcm = '{: 0.4F}'.format(cm / mxcm)
            for mkr in data[lg][cm]:
                curr_dat = data[lg][cm][mkr]
                try: rp1 = '{: 0.4F}'.format(curr_dat[4] / float(ftm_ln[curr_dat[3]])) 
                except: rp1 = 'NA'
                try: rp2 = '{: 0.4F}'.format(curr_dat[5] / float(ftm_ln[curr_dat[3]]))
                except: rp2 = 'NA'
                o.write('{0[0]}\t{0[1]}\t{0[2]}\t{0[3]}\t{0[4]}\t{0[5]}\t'.format(curr_dat))
                o.write(rcm + '\t' + rp1 + '\t' + rp2 + '\n')
    o.close()

def load_reflen():
    """
    Loads lengths of scaffolds into dict
    """
    dic = {}
    for r in reference.Reference(handle=reffa):
        scf = int(r.scaffold.replace('scaffold_',''))
        dic[scf] = r.len()
    return dic

if __name__=='__main__':
    start = time.time()
    if len(sys.argv) != 5:
        print main.__doc__
        sys.exit()
    inbest = sys.argv[1]
    linkin = sys.argv[2]
    reffa = sys.argv[3]
    outfil = sys.argv[4]
    main()
    print 'Finished; Took {: 0.4F} seconds.'.format(time.time()-start)
