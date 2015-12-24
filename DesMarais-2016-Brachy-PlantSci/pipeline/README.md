# Main Analysis Pipeline

## Step 00: Filter Reads

Raw SOLiD reads were filtered using the `ShortCSFilter.py` script available in this repository.

```
# Run filter command-line
ShortCSFilter.py <input.csfasta> <input.csqual> <filtered.fastq>
```

## Step 01: Align Reads

Filtered reads were aligned to the _B. distachyon_ genome v3.0 using [SHRiMP](http://compbio.cs.toronto.edu/shrimp/) (v2.2.3).
Alignments were filtered to remove unmapped reads and low quality (< 10) alignments using [sambamba](http://lomereiter.github.io/sambamba/)
(v0.5.8).

```
# Alignment and filtering command-line
gmapper-cs -N 4 -L <reference> --sam-header-rg <sam_header.txt> <filtered.fastq> \
    | sambamba view -S -o <filtered.bam> -f bam -F "not unmapped and mapping_quality > 10" /dev/stdin && \
    sambamba sort -o <sorted.bam> -t 4 -m 2GB --tmpdir=<tmp_directory> <filtered.bam>
```
