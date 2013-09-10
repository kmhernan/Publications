#!/usr/bin/env python
import sys
import reference

def main():
    """
    ---------------------------------------------------------------------------
    AUTHOR: Kyle Hernandez
    EMAIL: kmhernan@utexas.edu

    Create a karyotype config file. Of the form:
    # Chromosomes
    chr - ID LABEL START END COLOR
    # CYTOGENETIC BANDS
    band PARID LABEL LABEL START END COLOR
    ---------------------------------------------------------------------------
    USAGE: python make_karyotype.py infil ref.fa karyotype.txt links.txt 
    """
    sizes = load_ref()
    dat   = load_data(sizes)
    print_ktype(dat)

def load_ref():
    """
    Loads the reference fasta file into a dictionary
    """
    dic = {}
    ref_init = reference.Reference(handle=reffil)
    for r in ref_init:
        scf = int(r.scaffold.replace('scaffold_',''))
        dic[scf] = r.len()
    return dic

def load_data(sizes):
    """
    Loads R plot data into dict
    """
    dat = {'fm': {}, 'fh': {}}
    fm_scf = []
    skip = [11, 12, 14]
    o = open(outlk, 'wb')
    count = 0
    for line in open(infil, 'rU'):
        if line.startswith('LG'): continue
        else:
            cols = line.rstrip().split('\t')
            fhid = int(cols[0])
            fhst = float(cols[1])
            if fhid not in dat['fh']: dat['fh'][fhid] = {}
            dat['fh'][fhid][fhst] = ''
            try:
                fmid = int(cols[3])
                if fmid not in skip:
                    fm_scf.append(fmid)
                    count += 1 
                    print_links(cols, o, count)
            except:
                pass 
    o.close()

    for s in set(fm_scf):
        dat['fm'][s] = sizes[s]

    return dat
         
def print_links(cols, o, ct):
    """
    Creates dic to print link file
    """
    row1 = 'link' + str(ct) + ' FH{1[0]} {0:.0F} {0:.0F}'.format(float(cols[1]), cols)
    row2 = 'link' + str(ct) + ' FM{0[3]} {0[4]} {0[4]}'.format(cols)
    o.write(row1 + '\n' + row2 + '\n') 

def print_ktype(dat):
    """
    Prints to karyotype format.
    """
    #COLOR = {1:
    with open(outkt, 'wb') as o:
        for c in sorted(dat['fh']):
            o.write('chr - FH{0} {0} 0 {1:.0F} {2}\n'.format(c, max(dat['fh'][c].keys()) + 1, "blue"))
        for d in sorted(dat['fm']):
            o.write('chr - FM{0} {0} 0 {1} {2}\n'.format(d, dat['fm'][d], "green"))

if __name__ == '__main__':
    if len(sys.argv) != 5:
        print main.__doc__
        sys.exit()
    infil = sys.argv[1]
    reffil = sys.argv[2]
    outkt = sys.argv[3]
    outlk = sys.argv[4]
    main()
    print "Finished..."
