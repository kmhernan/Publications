# Main Analysis Pipeline

## Step 00: Filter Reads

Raw SOLiD reads were filtered using the `ShortCSFilter.py` script available in this repository.

```bash
# Run filter command-line
ShortCSFilter.py $input_fasta $input_csqual $filtered_fastq
```

## Step 01: Align Reads

Filtered reads were aligned to the _B. distachyon_ genome v3.0 using [SHRiMP](http://compbio.cs.toronto.edu/shrimp/) (v2.2.3).
Alignments were filtered to remove unmapped reads and low quality (< 10) alignments and sorted using [sambamba](http://lomereiter.github.io/sambamba/)
(v0.5.8). Finally, readgroup information was added to the bam files using [Picard](http://broadinstitute.github.io/picard/) (v1.138).

```bash
# Alignment and filtering command-line
gmapper-cs -N 4 -L $reference --sam-header-rg $sam_header_file $filtered_fastq \
    | sambamba view -S -o $filtered_bam -f bam -F "not unmapped and mapping_quality > 10" /dev/stdin \
    && sambamba sort -o $sorted_bam -t 4 -m 2GB --tmpdir=$tmp_directory $filtered.bam
    
# Add readgroups
java -Xmx1G -Xmx2G -XX:ParallelGCThreads=2 -jar $picard AddOrReplaceReadGroups I=$sorted_bam \
    O=$rg_bam \
    ID=${sample}.brachyrad \
    LB=${sample} \
    PL=solid \
    PU=${sample}.solid \
    SM=${sample}
    CN=UTGSF \
    CREATE_INDEX=true
```

## Step 02: GATK Processing

Filtered and sorted alignments were processed according to the 
[GATK's Best Practices](https://www.broadinstitute.org/gatk/guide/best-practices.php) (v3.4-46). 

```

```

Since there is currently
no source for high quality SNPs in the v3.0 genome, we followed 
