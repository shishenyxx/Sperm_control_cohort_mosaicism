# Snakemake pipeline for haplotypecaller for heterozygous SNVs and indels from the whole genome sequencing data

The pipeline is wirtten by Xinxu and Xiaoxu Yang

----------------------------

## Before start, make sure you have:
#### [java](https://www.java.com/en/download/help/linux_x64_install.xml) for Linux.
#### [GATK4](https://github.com/broadgsa/gatk/releases) version 4 for this pipeline.

----------------------------

## Input:
### Below are headers of the input file list
#### SAMPLE_NAME
User defined name for the sample.
#### Sex
The Sex of the sample 
#### BAM
The path to the input Dragen bam.
#### BAM_NAME
The ID of the sample from the bam file.
----------------------------

## Config files:
### Below are files you need to prepare for the annotation scripts, saved in the file snake_conf.yaml
#### n_intervals
number of intervals to split
#### input_paths
Path to the input file list.
#### out_dir
Path to the output directory.
#### scratch_dir
Path to the scratch files. Note that the number of temporary files will be euqal to two- or three-times the number of total listed variants.

#### java
Path to your java jre.
#### gatk4
Path to your GenomeAnalysisTK.jar (gatk4)

#### dbsnp
Snp list of your dbsnp file in vcf format (corresponding to your reference genome file).
#### ref_fasta
Your reference genome.
#### repetitive_region_bed
mills indel vcf file (corresponding to your reference genome file).
#### ref_bed
Length of each chromosome in the reference genome.



----------------------------

## Output:
The final output bam is in the output folder.


