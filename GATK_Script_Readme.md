# GATK Workflow Script for Alignment, Variant Calling and Quality Score Recalibration

## Overview
This script performs variant calling on paired-end whole-genome sequencing data using the GATK (Genome Analysis Toolkit). It aligns FASTQ files to a reference genome, processes the resulting BAM files, and applies variant quality score recalibration (VQSR) to produce high-quality variant calls in VCF format.

## Prerequisites

I am using conda environment, see how to create the identical one here:
[README.md](README.md)



### Resources

See the downloading instructions in [README.md](README.md)
- Reference genome: `Homo_sapiens_assembly38.fasta`
- Known variants for base recalibration and VQSR:
  - dbSNP: `Homo_sapiens_assembly38.dbsnp138.vcf`
  - Known indels: `Homo_sapiens_assembly38.known_indels.vcf.gz`
  - Mills and 1000G gold standard indels: `Mills_and_1000G_gold_standard.indels.hg38.vcf.gz`
  - HapMap: `hapmap_3.3.hg38.vcf.gz`
  - Omni: `1000G_omni2.5.hg38.vcf.gz`
  - 1000G Phase 1 high-confidence SNPs: `1000G_phase1.snps.high_confidence.hg38.vcf.gz`
- Target intervals: `wgs_calling_regions.hg38.interval_list`


## Script documentation
### Environment Setup

1. Load necessary modules:
   ```bash
   module purge
   module load gcc-libs
   ```
2. Activate the Conda environment:
   ```bash
   conda activate gatk_env
   ```
3. Ensure the script runs in the working directory:
   ```bash
   cd /home/skgtrk2/Scratch/<WORK_DIR>
   ```

## Workflow

### Step 1: Alignment with BWA-MEM
Align paired-end FASTQ files to the reference genome using BWA-MEM:
```bash
bwa mem -t 16 -R "@RG\tID:$SAMPLE\tSM:$SAMPLE\tPL:ILLUMINA" "$REFERENCE" "$FASTQ1" "$FASTQ2" | \
   samtools view -Sb - > "$OUTPUT_DIR/$SAMPLE.bam"
```

### Step 2: BAM Sorting
Sort the BAM file using Samtools:
```bash
samtools sort -@ 16 -o "$OUTPUT_DIR/$SAMPLE.sorted.bam" "$OUTPUT_DIR/$SAMPLE.bam"
rm "$OUTPUT_DIR/$SAMPLE.bam"
```

### Step 3: Mark Duplicates
Identify and mark duplicate reads:
```bash
gatk MarkDuplicates \
   -I "$OUTPUT_DIR/$SAMPLE.sorted.bam" \
   -O "$OUTPUT_DIR/$SAMPLE.marked_duplicates.bam" \
   -M "$OUTPUT_DIR/$SAMPLE.duplicate_metrics.txt" \
   --TMP_DIR "$OUTPUT_DIR/tmp"
```
Index the BAM file:
```bash
samtools index "$OUTPUT_DIR/$SAMPLE.marked_duplicates.bam"
```

### Step 4: Base Quality Score Recalibration (BQSR)
#### Generate Recalibration Table
```bash
gatk BaseRecalibrator \
   -I "$OUTPUT_DIR/$SAMPLE.marked_duplicates.bam" \
   -R "$REFERENCE" \
   --known-sites "$DBSNP" \
   --known-sites "$KNOWN_INDELS" \
   --known-sites "$MILLS_INDELS" \
   -O "$OUTPUT_DIR/$SAMPLE.recal_data.table" \
   --tmp-dir "$OUTPUT_DIR/tmp"
```
#### Apply Recalibration
```bash
gatk ApplyBQSR \
   -I "$OUTPUT_DIR/$SAMPLE.marked_duplicates.bam" \
   -R "$REFERENCE" \
   --bqsr-recal-file "$OUTPUT_DIR/$SAMPLE.recal_data.table" \
   -O "$OUTPUT_DIR/$SAMPLE.recalibrated.bam" \
   --tmp-dir "$OUTPUT_DIR/tmp"
```

### Step 5: Variant Calling
Generate a GVCF file:
```bash
gatk HaplotypeCaller \
   -R "$REFERENCE" \
   -I "$OUTPUT_DIR/$SAMPLE.recalibrated.bam" \
   -O "$OUTPUT_DIR/$SAMPLE.g.vcf.gz" \
   --dbsnp "$DBSNP" \
   -ERC GVCF \
   -L "$WGS_INTERVALS" \
   --native-pair-hmm-threads 16 \
   --tmp-dir "$OUTPUT_DIR/tmp"
```

### Step 6: Genotyping
```bash
gatk GenotypeGVCFs \
   -R "$REFERENCE" \
   -V "$OUTPUT_DIR/$SAMPLE.g.vcf.gz" \
   -O "$OUTPUT_DIR/$SAMPLE.vcf.gz" \
   --tmp-dir "$OUTPUT_DIR/tmp"
```

### Step 7: Variant Quality Score Recalibration (VQSR)
#### VQSR for SNPs
Generate recalibration file:
```bash
gatk VariantRecalibrator \
   -V "$OUTPUT_DIR/$SAMPLE.vcf.gz" \
   -R "$REFERENCE" \
   --resource:hapmap,known=false,training=true,truth=true,prior=15.0 "$HAPMAP" \
   --resource:omni,known=false,training=true,truth=true,prior=12.0 "$OMNI" \
   --resource:1000G,known=false,training=true,truth=false,prior=10.0 "$PHASE1_1000G" \
   --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 "$DBSNP" \
   -an DP -an QD -an FS -an SOR -an MQRankSum -an ReadPosRankSum \
   -mode SNP \
   -O "$OUTPUT_DIR/$SAMPLE.SNP.recal" \
   --tranches-file "$OUTPUT_DIR/$SAMPLE.SNP.tranches" \
   --rscript-file "$OUTPUT_DIR/$SAMPLE.SNP.plots.R" \
   --tmp-dir "$OUTPUT_DIR/tmp"
```
Apply recalibration:
```bash
gatk ApplyVQSR \
   -V "$OUTPUT_DIR/$SAMPLE.vcf.gz" \
   -O "$OUTPUT_DIR/$SAMPLE.SNP.vcf.gz" \
   --recal-file "$OUTPUT_DIR/$SAMPLE.SNP.recal" \
   --tranches-file "$OUTPUT_DIR/$SAMPLE.SNP.tranches" \
   --truth-sensitivity-filter-level 99.5 \
   -mode SNP \
   --tmp-dir "$OUTPUT_DIR/tmp"
```
#### VQSR for Indels
Generate recalibration file:
```bash
gatk VariantRecalibrator \
   -V "$OUTPUT_DIR/$SAMPLE.SNP.vcf.gz" \
   -R "$REFERENCE" \
   --resource:mills,known=true,training=true,truth=true,prior=12.0 "$MILLS_INDELS" \
   --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 "$DBSNP" \
   -an DP -an FS -an MQRankSum -an ReadPosRankSum -an SOR \
   -mode INDEL \
   -O "$OUTPUT_DIR/$SAMPLE.INDEL.recal" \
   --tranches-file "$OUTPUT_DIR/$SAMPLE.INDEL.tranches" \
   --rscript-file "$OUTPUT_DIR/$SAMPLE.INDEL.plots.R" \
   --tmp-dir "$OUTPUT_DIR/tmp"
```
Apply recalibration:
```bash
gatk ApplyVQSR \
   -V "$OUTPUT_DIR/$SAMPLE.SNP.vcf.gz" \
   -O "$OUTPUT_DIR/$SAMPLE.SNP.INDEL.vcf.gz" \
   --recal-file "$OUTPUT_DIR/$SAMPLE.INDEL.recal" \
   --tranches-file "$OUTPUT_DIR/$SAMPLE.INDEL.tranches" \
   --truth-sensitivity-filter-level 99.0 \
   -mode INDEL \
   --tmp-dir "$OUTPUT_DIR/tmp"
```

