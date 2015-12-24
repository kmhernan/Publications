#!/usr/bin/env python
"""
Kyle Hernandez
khernandez@bsd.uchicago.edu
Center for Research Informatics
University of Chicago
"""
import sys

hpoly = 0.2
minq = 20.0
bases = ['0', '1', '2', '3']

def usage():
    """Prints usage to the screen"""
    print """
-------------------------------------------------------------------------------
Author: Kyle Hernandez <khernandez@bsd.uchicago.edu>

Description: Filters SOLiD CS fasta/qual files, converts them to fastq files, 
and writes out filtering statistics to a file. The statistics file will be the
same as the <filtered.fq> file with an '.trim_stats' added to the end.

Filters: There are hard-coded filters applied in this order:
	- Reads with missing base calls ("missing")
	- Reads where the quality length > read length ("lengths")
	- Reads with a mean quality score < 20.0 ("lowqual")
	- Reads with homopolymer lengths > 0.2 * read length ("homopolymer")  

Usage:
	ShortCSFilter.py <raw.fa> <raw.qual> <filtered.fq>
-------------------------------------------------------------------------------
"""

def process_files(in_fasta, in_qual, o_fq):
    """
    Performs filtering, counting, and conversion
    @param in_fasta - input raw csfasta file
    @param in_qual - input raw csqual file
    @param o_fq - output filtered fastq file
    @returns dictionary of filtering statistics
    """
    q_fh = open(in_qual, 'rU')
    f_fh = open(in_fasta, 'rU')
    o    = open(o_fq, 'wb')
    cts  = {'total': 0, 'lengths': 0, 'homopolymer': 0, 'lowqual': 0, 'missing': 0, 'passed': 0} 

    # Loop over both csfasta and csqual files at the same time
    try:
        while True:
            # csfasta info
            f_head = f_fh.readline().rstrip()
            if not f_head: break
            f_seq  = f_fh.readline().rstrip()
   
            # csqual info 
            q_head  = q_fh.readline().rstrip()
            q_qual  = q_fh.readline().rstrip()

            # Get the quality array
            q_array = map(int, q_qual.split(' '))
  
            # Total sequence records 
            cts['total'] += 1 

            # Assert the csfasta and csqual identifiers are the same
            assert f_head == q_head, 'ERROR!! Line {0}. Indentifier mismatch:\nfasta: {1}\ncsqual: {2}'.format(
                cts['total'], f_head, q_head)

            # Missing filter
            if '.' in f_seq: 
                cts['missing'] += 1
                continue
        
            # Lengths filter 
            if len(q_array) > len(f_seq) - 1:
                cts['lengths'] += 1
                continue 
        
            # Get average quality and filter
            mean_q   = sum(q_array) / float(len(q_array)) 
            if mean_q < minq:
                cts['lowqual'] += 1
                continue

            # Trim reads to be of the same length as qual 
            fix_seq  = f_seq[:len(q_array) + 1]

            # Filter homopolymers
            for b in bases:
                curr = b*int(len(q_array) * hpoly) 
                if curr in fix_seq:
                    cts['homopolymer'] += 1
                    continue

            # Passed filters
            cts['passed'] += 1

            # Convert to fastq and write
            o.write('@' + f_head[1:] + '\n')
            o.write(fix_seq + '\n')
            o.write('+\n')
            ascii_qual = [chr(i+33) if i <=40 else 'I' for i in q_array]
            o.write(''.join(ascii_qual) + '\n')
    finally:
        q_fh.close()
        f_fh.close()
        o.close()

    # Close handles
    q_fh.close()
    f_fh.close()
    o.close()

    # Return stats
    return cts

if __name__ == '__main__':
    if len(sys.argv) != 4:
        usage() 
        sys.exit(1)

    # Input files
    in_fasta = sys.argv[1]
    in_qual  = sys.argv[2]
    o_fq     = sys.argv[3]

    # Filter, write fastq, and get stats
    cmap = process_files(in_fasta, in_qual, o_fq)

    # Write stats
    o_stats  = o_fq + '.trim_stats'
    keys = ['total', 'lengths', 'missing', 'lowqual', 'homopolymer', 'passed']
    with open(o_stats, 'wb') as o:
        for k in keys:
            o.write(k + '\t' + str(cmap[k]) + '\n') 
