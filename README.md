# GATK_Myriad

The actual script is runGATK.sh

### Preparation instructions
I created a separate conda environment where I put GATK and samtools and a few other dependencies


### DOwnload all GATK dependencies from google clour Broad Institute
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
```
