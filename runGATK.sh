#!/bin/bash -l
#$ -l h_rt=24:00:00
#$ -l mem=35G
#$ -l tmpfs=35G
#$ -N GATK
#$ -pe smp 16

module purge
module load gcc-libs

conda activate gatk_env
cd /home/skgtrk2/Scratch/ICGNMD_TEMP

# Define paths to FASTQ files
FASTQ1="/home/skgtrk2/Scratch/ICGNMD_TEMP/IC_KUC_00124_WGS_macrogen_R1.fastq.gz"
FASTQ2="/home/skgtrk2/Scratch/ICGNMD_TEMP/IC_KUC_00124_WGS_macrogen_R2.fastq.gz"
OUTPUT_DIR="gatk_output"
SAMPLE="IC_KUC_00124"

# Define paths to resources in your bundle directory
RESOURCES_DIR="/home/skgtrk2/Scratch/GATK_resource/v0"
REFERENCE="$RESOURCES_DIR/Homo_sapiens_assembly38.fasta"
DBSNP="$RESOURCES_DIR/Homo_sapiens_assembly38.dbsnp138.vcf"
KNOWN_INDELS="$RESOURCES_DIR/Homo_sapiens_assembly38.known_indels.vcf.gz"
MILLS_INDELS="$RESOURCES_DIR/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
HAPMAP="$RESOURCES_DIR/hapmap_3.3.hg38.vcf.gz"
OMNI="$RESOURCES_DIR/1000G_omni2.5.hg38.vcf.gz"
PHASE1_1000G="$RESOURCES_DIR/1000G_phase1.snps.high_confidence.hg38.vcf.gz"
WGS_INTERVALS="$RESOURCES_DIR/wgs_calling_regions.hg38.interval_list"


# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"
### STEP 1: Align FASTQ files to reference genome with BWA-MEM and convert to BAM ###
echo "starting bwa mem"
bwa mem -t 16 -R "@RG\tID:$SAMPLE\tSM:$SAMPLE\tPL:ILLUMINA" "$REFERENCE" "$FASTQ1" "$FASTQ2" | \
   samtools view -Sb - > "$OUTPUT_DIR/$SAMPLE.bam"
echo "done bwa mem, staring samtools sorting"
# Sort BAM file
samtools sort -@ 16 -o "$OUTPUT_DIR/$SAMPLE.sorted.bam" "$OUTPUT_DIR/$SAMPLE.bam"
rm "$OUTPUT_DIR/$SAMPLE.bam"
echo "done samtools sorting, starting marking duplicates"
### STEP 2: Mark duplicates ###

gatk MarkDuplicates \
   -I "$OUTPUT_DIR/$SAMPLE.sorted.bam" \
   -O "$OUTPUT_DIR/$SAMPLE.marked_duplicates.bam" \
   -M "$OUTPUT_DIR/$SAMPLE.duplicate_metrics.txt" \
   --TMP_DIR "$OUTPUT_DIR/tmp"

# Index the BAM
samtools index "$OUTPUT_DIR/$SAMPLE.marked_duplicates.bam"
echo "staring base recalibrator 1"
### STEP 3: Base Quality Score Recalibration (BQSR) ###
# Step 3a: Generate recalibration table
gatk BaseRecalibrator \
   -I "$OUTPUT_DIR/$SAMPLE.marked_duplicates.bam" \
   -R "$REFERENCE" \
   --known-sites "$DBSNP" \
   --known-sites "$KNOWN_INDELS" \
   --known-sites "$MILLS_INDELS" \
   -O "$OUTPUT_DIR/$SAMPLE.recal_data.table" \
   --tmp-dir "$OUTPUT_DIR/tmp" \
   --verbosity INFO

echo "apply it"
# Step 3b: Apply BQSR
gatk ApplyBQSR \
   -I "$OUTPUT_DIR/$SAMPLE.marked_duplicates.bam" \
   -R "$REFERENCE" \
   --bqsr-recal-file "$OUTPUT_DIR/$SAMPLE.recal_data.table" \
   -O "$OUTPUT_DIR/$SAMPLE.recalibrated.bam" \
   --tmp-dir "$OUTPUT_DIR/tmp"


echo "call haplotype caller"
### STEP 4: Variant Calling with HaplotypeCaller in GVCF mode ###
gatk HaplotypeCaller \
   -R "$REFERENCE" \
   -I "$OUTPUT_DIR/$SAMPLE.recalibrated.bam" \
   -O "$OUTPUT_DIR/$SAMPLE.g.vcf.gz" \
   --dbsnp "$DBSNP" \
   -ERC GVCF \
   -L "$WGS_INTERVALS" \
   --native-pair-hmm-threads 16 \
   --tmp-dir "$OUTPUT_DIR/tmp"

echo "done haplotype caller, start genotype"
### STEP 5: Joint Genotyping (if you have multiple samples, include all GVCFs here) ###
gatk GenotypeGVCFs \
   -R "$REFERENCE" \
   -V "$OUTPUT_DIR/$SAMPLE.g.vcf.gz" \
   -O "$OUTPUT_DIR/$SAMPLE.vcf.gz" \
   --tmp-dir "$OUTPUT_DIR/tmp"

echo "recalbrate for snps"
### STEP 6: Apply VQSR (Variant Quality Score Recalibration) ###
# VQSR for SNPs
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

echo "apply it"
gatk ApplyVQSR \
   -V "$OUTPUT_DIR/$SAMPLE.vcf.gz" \
   -O "$OUTPUT_DIR/$SAMPLE.SNP.vcf.gz" \
   --recal-file "$OUTPUT_DIR/$SAMPLE.SNP.recal" \
   --tranches-file "$OUTPUT_DIR/$SAMPLE.SNP.tranches" \
   --truth-sensitivity-filter-level 99.5 \
   -mode SNP \
   --tmp-dir "$OUTPUT_DIR/tmp"

echo "recalibrator for indels"
# VQSR for INDELs
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
echo "apply it"
gatk ApplyVQSR \
   -V "$OUTPUT_DIR/$SAMPLE.SNP.vcf.gz" \
   -O "$OUTPUT_DIR/$SAMPLE.SNP.INDEL.vcf.gz" \
   --recal-file "$OUTPUT_DIR/$SAMPLE.INDEL.recal" \
   --tranches-file "$OUTPUT_DIR/$SAMPLE.INDEL.tranches" \
   --truth-sensitivity-filter-level 99.0 \
   -mode INDEL \
   --tmp-dir "$OUTPUT_DIR/tmp"
