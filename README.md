# Sperm control cohort mosaicism
This repository collects pipelines, codes, and some intermediate results for the study of mosaic SNV/Indels for sperm, blood, and other samples of a control cohort. Raw data of this study is available on SRA under [PRJNA660493](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA660493/) and [PRJNA588332](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA588332/).

### 1. Pipelines for the process of whole-genome sequencing data
#### 1.1 Pipelines for WGS data process and quality control
#### 1.2 Codes for the population origin analysis
#### 1.3 Pipelines for mosaic SNV/indel calling and variant annotations
[Pipelines](https://github.com/shishenyxx/Adult_brain_somatic_mosaicism/tree/master/pipelines/WGS_SNV_indel_calling_pipeline/Mutect2_PM_Strelka2) for MuTect2 (paired mode) and Strelka2 (somatic mode) variant calling from WGS data

Pipelines for [MuTect2 (single mode)](https://github.com/shishenyxx/Adult_brain_somatic_mosaicism/tree/master/pipelines/WGS_SNV_indel_calling_pipeline/Mutect2_single_mode), followed by [MosaicForecast](https://github.com/shishenyxx/Adult_brain_somatic_mosaicism/tree/master/pipelines/WGS_SNV_indel_calling_pipeline/MosaicForecast_pipeline), and the [variant annotation pipeline](https://github.com/shishenyxx/PASM/tree/master/Snakemake_pipeline).

PBS script for [MosaicHunter (single mode)](https://github.com/shishenyxx/Adult_brain_somatic_mosaicism/tree/master/pipelines/WGS_SNV_indel_calling_pipeline/MosaicHunter_single_mode_pipeline), followed by the [variant annotation pipeline](https://github.com/shishenyxx/PASM/tree/master/Snakemake_pipeline).

After variant calling from different strategies, variants were annotated and filtered by [a python script](https://github.com/shishenyxx/Adult_brain_somatic_mosaicism/blob/master/pipelines/WGS_SNV_indel_calling_pipeline/Filter_and_annotate_candidate_mosaic_variants.py) and positive mosaic variants as well as the corresponding tissue and additional information were annotated.

### 2. Pipelines for the process of Targeted Amplicon Sequencing (TAS)
#### 2.1 Pipelines for TAS data alignment and processing
#### 2.2 Pipelines for AF quantification and variant annotations

### 3. Pipelines for the data analysis, variant filtering, comprehensive annotations, and statistical analysis
#### 3.1 Pipelines for mosaic variant determination, annotations, and plotting
#### 3.2 Pipelines for statistically analysis, QC, and the related plotting
#### 3.3 Codes for statistical modeling
[Codes](https://github.com/shishenyxx/Sperm_control_cohort_mosaicism/blob/master/Mutation_accumulation_model/fit_afs.py) for the estimation of accumulation of mutations through a stepwise exponential regression regression model.
