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

```bash
# Get targets
java -Xmx2G -Xmx4G -XX:ParallelGCThreads=2 -jar $gatk -T RealignerTargetCreator \
    -R $reference \
    -I $rg_bam \
    -nt 2 \
    -o $targets_intervals

# Realign
java -Xmx2G -Xmx4G -XX:ParallelGCThreads=2 -jar $gatk -T IndelRealigner \
    -R $reference \
    -I $rg_bam \
    -o $realigned_bam \
    --targetIntervals $targets_intervals
```

Since there is currently no source for high quality SNPs in the v3.0 genome, we first ran all of the indel realigned bams through 
a single run of the GATK's UnifiedGenotyper. Then, we pulled out the top quality SNPs and used those to run the base quality
recalibration step.

```bash
# Run UnifiedGenotyper
java -Xms10G -Xmx25G -XX:ParallelGCThreads=2 -jar $gatk -T UnifiedGenotyper \
    -nct 4 -nt 4 \
    -I $bam_file_list \
    -stand_call_conf 20.0 \
    -stand_emit_conf 20.0 \
    -glm SNP \
    -gt_mode DISCOVERY \
    -R $ref -L $chrom -o $raw_vcf

# Run basic filters
java -Xmx8G -Xmx10G -XX:ParallelGCThreads=2 -jar $gatk -T VariantFiltration \
    -V $raw_vcf \
    -R $reference \
    -L $chrom \
    -o $filtered_vcf \
    --filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
    --filterName "BigFilter"

# Get top quality SNPs using the GetHighQualSNPs.py script in this repository
GetHighQualSNPs.py $filtered_vcf $top_vcf

# Recalibration
java -Xms2G -Xmx4G -jar $gatk -T BaseRecalibrator \
    -R $reference \
    -I $realigned_bam \
    -knownSites $top_vcf
    -o $recal_table

# Apply Recalibration
java -Xms2G -Xmx4G -jar $gatk -T PrintReads \
    -R $reference \
    -BQSR $recal_table \
    -I $realigned_bam \
    -o $recalibrated.bam
```
