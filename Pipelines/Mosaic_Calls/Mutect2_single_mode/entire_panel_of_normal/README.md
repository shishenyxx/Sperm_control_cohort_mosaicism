# Snakemake pipeline for the single mode of MuTect2 using a "entire" panel of normal strategy

The pipeline is wirtten by Xin Xu, with help form Xiaoxu Yang

----------------------------

## Before start, make sure you have:
#### [java](https://www.java.com/en/download/help/linux_x64_install.xml) for Linux.
#### [GATK4](https://github.com/broadgsa/gatk/releases) version 4 for this pipeline.

----------------------------

## Input:
### Below are headers of the input file list
#### sample	
Unique IDs for any input bam file.
#### tumor	
The ID of the sample from the bam file.
#### tumor_path
The path to the input bam.

----------------------------

## Config files:
### Below are files you need to prepare for the pipeline, saved in the file snake_conf.yaml
#### n_intervals
number of intervals to split, for cluster job submission
#### input_files
Path to the input file list for the panel of normals.
#### out_dir
Path to the output directory.
#### scratch_dir
Path to the scratch files. Note that the number of temporary files will be euqal to two- or three-times the number of total listed variants.

#### java
Path to your java jre.
#### gatk4
Path to your GenomeAnalysisTK.jar (GATK4)

#### ref_fasta
Your reference genome.
#### bed_file 
Length for each chromosome (corresponding to your reference genome file).
#### gnomad
Path to the gnomAD allelic frequencies in vcf,gz format (corresponding to your reference genome file).

----------------------------

## Output:
The final output bam is in the results folder.


