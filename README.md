# GATK_Myriad

The actual script is runGATK.sh

The script and all its steps are documented in detail [GATK_Script_Readme.md](GATK_Script_Readme.md)

### Preparation instructions
I created a separate conda environment where I put GATK, samtools and a bwa
```
conda create -n gatk_env
conda activate gatk_env
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda install gatk4
conda install -c bioconda samtools bwa
conda install -c r r-base r-ggplot2

#to check what is installed in the activated env
conda env export --from-history
```

### DOwnload all GATK dependencies from google cloud Broad Institute
```
conda create -n gcloud_env python=3.8 -y
conda activate gcloud_env
conda install -c conda-forge google-cloud-sdk
conda install -c conda-forge crcmod
#check that crcmod is used
python -c "import crcmod; print(crcmod._usingExtension)"
#it was not install properly, trying with pip
pip uninstall -y crcmod
pip install --no-binary :all: crcmod
python -c "import crcmod; print(getattr(crcmod, '_usingExtension', 'Attribute not found'))"
#worked fine with pip

#this works via a browser on another machine, follow on screen instructions
gcloud auth login --no-launch-browser
gsutil ls gs://genomics-public-data/resources/broad/hg38/v0
gsutil -m cp -r gs://genomics-public-data/resources/broad/hg38/v0 .
#to skip existing files
gsutil -m cp -r -n gs://genomics-public-data/resources/broad/hg38/v0 .
#this crashed a few times, so trying various flags
gsutil -D cp -r -n gs://genomics-public-data/resources/broad/hg38/v0 .
```
This will result in these files being downloaded:
```
cd /home/skgtrk2/Scratch/GATK_resource/v0
ls -lth

total 33G
-rw-r--r--  1 skgtrk2 skgts1  11G Nov  4 18:09 Homo_sapiens_assembly38.dbsnp138.vcf
-rw-r--r--  1 skgtrk2 skgts1  12G Nov  4 17:53 1000G.phase3.integrated.sites_only.no_MATCHED_REV.hg38.vcf
-rw-r--r--  1 skgtrk2 skgts1 3.1G Nov  4 17:50 Homo_sapiens_assembly38.fasta
-rw-r--r--  1 skgtrk2 skgts1 3.0G Nov  4 17:50 Homo_sapiens_assembly38.fasta.64.bwt
-rw-r--r--  1 skgtrk2 skgts1 1.8G Nov  4 17:49 1000G_phase1.snps.high_confidence.hg38.vcf.gz
-rw-r--r--  1 skgtrk2 skgts1 1.5G Nov  4 17:49 Homo_sapiens_assembly38.fasta.64.sa
-rw-r--r--  1 skgtrk2 skgts1 768M Nov  4 17:49 Homo_sapiens_assembly38.fasta.64.pac
-rw-r--r--  1 skgtrk2 skgts1  60M Nov  4 17:49 hapmap_3.3.hg38.vcf.gz
-rw-r--r--  1 skgtrk2 skgts1  51M Nov  4 17:49 1000G_omni2.5.hg38.vcf.gz
-rw-r--r--  1 skgtrk2 skgts1  59M Nov  4 17:49 Homo_sapiens_assembly38.known_indels.vcf.gz
-rw-r--r--  1 skgtrk2 skgts1  12M Nov  4 17:48 Homo_sapiens_assembly38.dbsnp138.vcf.idx
-rw-r--r--  1 skgtrk2 skgts1  20M Nov  4 17:48 Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
-rw-r--r--  1 skgtrk2 skgts1 1.5M Nov  4 17:48 1000G_omni2.5.hg38.vcf.gz.tbi
-rw-r--r--  1 skgtrk2 skgts1 2.1M Nov  4 17:48 1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi
-rw-r--r--  1 skgtrk2 skgts1 9.3M Nov  4 17:48 1000G.phase3.integrated.sites_only.no_MATCHED_REV.hg38.vcf.idx
-rw-r--r--  1 skgtrk2 skgts1 3.0M Nov  4 17:48 Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz
-rw-r--r--  1 skgtrk2 skgts1 412K Nov  4 17:48 Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz.tbi
-rw-r--r--  1 skgtrk2 skgts1 1.5M Nov  4 17:48 hapmap_3.3.hg38.vcf.gz.tbi
-rw-r--r--  1 skgtrk2 skgts1 569K Nov  4 17:48 Homo_sapiens_assembly38.dict
-rw-r--r--  1 skgtrk2 skgts1 477K Nov  4 17:48 Homo_sapiens_assembly38.fasta.64.alt
-rw-r--r--  1 skgtrk2 skgts1  20K Nov  4 17:48 Homo_sapiens_assembly38.fasta.64.amb
-rw-r--r--  1 skgtrk2 skgts1 445K Nov  4 17:48 Homo_sapiens_assembly38.fasta.64.ann
-rw-r--r--  1 skgtrk2 skgts1 158K Nov  4 17:48 Homo_sapiens_assembly38.fasta.fai
-rw-r--r--  1 skgtrk2 skgts1 1.5M Nov  4 17:48 Homo_sapiens_assembly38.known_indels.vcf.gz.tbi
-rw-r--r--  1 skgtrk2 skgts1 1.5M Nov  4 17:48 Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi
-rw-r--r--  1 skgtrk2 skgts1 586K Nov  4 17:48 wgs_calling_regions.hg38.interval_list
drwxr-xr-x 52 skgtrk2 skgts1 4.0K Nov  4 17:28 scattered_calling_intervals

```
