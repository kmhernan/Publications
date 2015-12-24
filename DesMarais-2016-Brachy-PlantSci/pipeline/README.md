# Main Analysis Pipeline

## Step 00: Filter Reads

Raw SOLiD reads were filtered using the `ShortCSFilter.py` script available in this repository.

```
# Run filter command-line skeleton
ShortCSFilter.py <input.csfasta> <input.csqual> <filtered.fastq>
```

## Step 01: Align Reads

Filtered reads were aligned to the _B. distachyon_ genome v3.0 using [SHRiMP](http://compbio.cs.toronto.edu/shrimp/) (v2.2.3).
Alignments were filtered to remove unmapped reads and low quality (< 10) alignments using [sambamba](http://lomereiter.github.io/sambamba/)
(v0.5.8).
