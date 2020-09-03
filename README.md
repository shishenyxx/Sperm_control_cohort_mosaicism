# Sperm control cohort mosaicism

This repository collects pipelines, codes, and some intermediate results for the study of mosaic SNV/Indels for sperm, blood, and other samples of a control cohort. Raw data of this study is available on SRA under [PRJNA660493](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA660493/) and [PRJNA588332](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA588332/).

### 1. Pipelines for the process of whole-genome sequencing data

#### 1.1 Pipelines for WGS data process and quality control

[Pipelines](https://github.com/shishenyxx/Sperm_control_cohort_mosaicism/tree/master/Pipelines/Preprocessing) for pre-processing of the bams.

Codes for [depth of coverage](https://github.com/shishenyxx/Sperm_control_cohort_mosaicism/blob/master/Plot/QC/Depth_of_coverage.r) and [insertsize distribution](https://github.com/shishenyxx/Sperm_control_cohort_mosaicism/blob/master/Plot/QC/Insert_size.r).

#### 1.2 Codes for the population origin analysis

[Pipeline](https://github.com/shishenyxx/Sperm_control_cohort_mosaicism/blob/master/Plot/Population_analysis/Whole_genome_variant_extraction_and_PCA.sh) for population analysis, and [codes](https://github.com/shishenyxx/Sperm_control_cohort_mosaicism/blob/master/Plot/Population_analysis/Plot_PCA.r) for plot.

#### 1.3 Pipelines for mosaic SNV/indel calling and variant annotations

[Pipelines](https://github.com/shishenyxx/Adult_brain_somatic_mosaicism/tree/master/pipelines/WGS_SNV_indel_calling_pipeline/Mutect2_PM_Strelka2) for MuTect2 (paired mode) and Strelka2 (somatic mode) variant calling from WGS data

Pipelines for MuTect2 (single mode) has a ["Leave One Out"](https://github.com/shishenyxx/Sperm_control_cohort_mosaicism/tree/master/Pipelines/Mosaic_Calls/Mutect2_single_mode/leave_one_out) version for the YA cohort, and a ["Full Panel of Normal"](https://github.com/shishenyxx/Adult_brain_somatic_mosaicism/tree/master/pipelines/WGS_SNV_indel_calling_pipeline/Mutect2_single_mode) version for the AA and ASD cohort. The MuTect2 (single mode) result is followed by [MosaicForecast](https://github.com/shishenyxx/Adult_brain_somatic_mosaicism/tree/master/pipelines/WGS_SNV_indel_calling_pipeline/MosaicForecast_pipeline), and the [variant annotation pipeline](https://github.com/shishenyxx/PASM/tree/master/Snakemake_pipeline).

[Codes] to plot the calls of different methods on simulated variants.

### 2. Pipelines for the process of Targeted Amplicon Sequencing (TAS)
#### 2.1 Pipelines for TAS data alignment and processing

[Pipelines](https://github.com/shishenyxx/Adult_brain_somatic_mosaicism/tree/master/pipelines/MPAS_and_snMPAS_processing_pipeline) for alignment, processing, and germline variant calling of TAS reads.

#### 2.2 Pipelines for AF quantification and variant annotations

[Pipelines](https://github.com/shishenyxx/PASM/tree/master/Snakemake_pipeline) for AF quantification and variant anntations.

[Codes](https://github.com/shishenyxx/Sperm_control_cohort_mosaicism/blob/master/Pipelines/cc_validation_table_sperm_abc_saliva_blood.py) to filter and annotate on TAS data.

### 3. Pipelines for the data analysis, variant filtering, comprehensive annotations, and statistical analysis
#### 3.1 Pipelines for mosaic variant determination, annotations, and plotting

After variant calling from different strategies, variants were annotated and filtered by [a python script](https://github.com/shishenyxx/Sperm_control_cohort_mosaicism/blob/master/Pipelines/Mosaic_Calls/control_cohort_complete_from_variant_table_MSMF03.py) and positive mosaic variants as well as the corresponding samples and additional information were annotated.

[Codes](https://github.com/shishenyxx/Adult_brain_somatic_mosaicism/tree/master/permutation) for permutation analysis from gnomAD and [codes](https://github.com/shishenyxx/Sperm_control_cohort_mosaicism/tree/master/Plot/cc_annotation_counts02.py) for plotting the permutation result.

UpSet plot is generated from an [online tool](https://vcg.github.io/upset/).

#### 3.2 Pipelines for statistically analysis, and the related plotting

[Codes](https://github.com/shishenyxx/Sperm_control_cohort_mosaicism/blob/master/Mutation_accumulation_model/fit_afs.py) for the estimation of accumulation of mutations through a stepwise exponential regression regression model.

[Codes]() for the analysis of accuracy of number of variants in different groups. 
