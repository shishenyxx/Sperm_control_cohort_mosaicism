# Snakemake pipeline for Base Quality Score Recalibration for the input bam. The pipeline is originally wirtten by Xinxu, and modified by Xiaoxu Yang, with help form Martin Breuss.

----------------------------

## Before start, make sure you have:
#### [java](https://www.java.com/en/download/help/linux_x64_install.xml) for Linux.
#### [GATK4](https://github.com/broadgsa/gatk/releases) GATK versioni 4 for this pipeline.

----------------------------

## Input:
### Below are headers of the input file list
#### SAMPLE_NAME
Unique IDs for any input bam file.
#### BAM_NAME
The path to the input bam.

----------------------------

## Config files:
### Below are files you need to prepare for the annotation scripts, saved in the file snake_conf.yaml
#### bam_paths
Path to the input file list.
#### out_dir
Path to the output directory.
#### scratch_dir
Path to the scratch files. Note that the number of temporary files will be euqal to two- or three-times the number of total listed variants.
#### ref_fasta
Your reference genome.

#### gatk4
Path to run your (java -jar GenomeAnalysisTK.jar)

#### mills_indel
vmills indel vcf file (corresponding to your reference genome file).
#### indels_1000g
List for common indels from the 1000 Genome Project in vcf format (corresponding to your reference genome file).
#### dbsnp
Snp list of your dbsnp file in vcf format (corresponding to your reference genome file).


----------------------------

## Output:
The final output bam is in the recaled_bams folder.


