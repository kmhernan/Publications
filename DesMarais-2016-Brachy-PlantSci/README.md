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

## Step 03: Detect SNPs with Freebayes

We used [Freebayes](https://github.com/ekg/freebayes) (v0.9.21-19-gc003c1e) to detect variants in the mapping population. All
multi-allelelic loci were removed and variants were normalized using [Vt](http://genome.sph.umich.edu/wiki/Vt) (v0.5).

```bash
# Run freebayes
freebayes -f $reference -L $recal_bam_list_file -v $raw_vcf \
    -m 10 -q 5 --genotype-qualities --use-mapping-quality \
    --site-selection-max-iterations 3 --genotyping-max-iterations 25 \
    --min-alternate-count 2 --min-alternate-qsum 40 \
    --genotype-variant-threshold 4

# Get biallelic loci and remove loci with QUAL < 15
cat $raw_vcf | awk -F"\t" '{if($1~/^#/){print $0}else if($4!~/,/ && $5!~/,/ && $6>=15){print $0}}' > $biallelic_vcf

# Normalize
vt normalize -o $biallelic_normalized_vcf -r $reference $biallelic_vcf
```

## Step 04: Marker Calling and Filtering

The biallelic, normalized variants from __Step 02__ were used to determine which parent contributed the alleles. Genotypes 
were converted to A/B/H calls and then filtered to remove likely uninformative markers.

```bash
# Call markers
FreebayesToMarkers.py $biallelic_normalized_vcf $raw_markers_txt

# Filter markers
FilterMarkers.py $raw_markers_txt $filtered_markers_txt
```

## Step 05: Updating Cui Marker Locations

Since the Cui markers were aligned to the older version of the _B. distachyon_ genome, the flanking sequence from the markers
were aligned to the new genome using [blastn](http://www.ncbi.nlm.nih.gov/books/NBK279690/) (v2.2.31). Then, the top unique hits
were extracted into a bed format for further analysis.

```bash
# Blast Cui markers to new assembly
$blastn -db $dbname -query $cui_markers_fasta \
-out $cui_markers_update_txt \
-outfmt 7 -num_alignments 10 -num_threads 4

# Get top unique hits
ParseBlast.py $cui_markers_update_txt $cui_markers_best_bed
```
