# Sperm control cohort mosaicism

This repository collects pipelines, codes, and some intermediate results for the study of mosaic SNV/Indels for sperm, blood, and other samples of a control cohort. Raw data of this study is available [here](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA660493/) and [here](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA588332/). The [first](https://trace.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA660493&o=acc_s%3Aa) and [second](https://trace.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA588332&o=acc_s%3Aa) dataset can be accessed through [SRA Run Selector](https://trace.ncbi.nlm.nih.gov/Traces/study/?).

<img src="https://user-images.githubusercontent.com/17311837/223877525-ce835f79-274a-432b-9bf2-3e9405c41499.png" alt="Sperm_Mosaic_Cover" width=50%> 

-----------------------------------

### 1. Pipelines for the process of whole-genome sequencing data

#### 1.1 Pipelines for WGS data process and quality control

[Pipelines](https://github.com/shishenyxx/Sperm_control_cohort_mosaicism/tree/master/Pipelines/Preprocessing) for pre-processing of the bams.

Codes for [depth of coverage](https://github.com/shishenyxx/Sperm_control_cohort_mosaicism/blob/master/Plot/QC/Depth_of_coverage.r) and [insertsize distribution](https://github.com/shishenyxx/Sperm_control_cohort_mosaicism/blob/master/Plot/QC/Insert_size.r).

#### 1.2 Codes for the population origin analysis

[Pipeline](https://github.com/shishenyxx/Sperm_control_cohort_mosaicism/blob/master/Plot/Population_analysis/Whole_genome_variant_extraction_and_PCA.sh) for population analysis, and [codes](https://github.com/shishenyxx/Sperm_control_cohort_mosaicism/blob/master/Plot/Population_analysis/Plot_PCA.r) for plot.

#### 1.3 Pipelines for mosaic SNV/indel calling and variant annotations

[Pipelines](https://github.com/shishenyxx/Adult_brain_somatic_mosaicism/tree/master/pipelines/WGS_SNV_indel_calling_pipeline/Mutect2_PM_Strelka2) for MuTect2 (paired mode) and Strelka2 (somatic mode) variant calling from WGS data

Pipelines for MuTect2 (single mode) has a ["Leave One Out"](https://github.com/shishenyxx/Sperm_control_cohort_mosaicism/tree/master/Pipelines/Mosaic_Calls/Mutect2_single_mode/leave_one_out) version for the YA cohort, and a ["Full Panel of Normal"](https://github.com/shishenyxx/Adult_brain_somatic_mosaicism/tree/master/pipelines/WGS_SNV_indel_calling_pipeline/Mutect2_single_mode) version for the AA and ASD cohort. The MuTect2 (single mode) result is followed by [MosaicForecast](https://github.com/shishenyxx/Adult_brain_somatic_mosaicism/tree/master/pipelines/WGS_SNV_indel_calling_pipeline/MosaicForecast_pipeline), and the [variant annotation pipeline](https://github.com/shishenyxx/PASM/tree/master/Snakemake_pipeline).

[Codes](https://github.com/shishenyxx/Sperm_control_cohort_mosaicism/blob/master/Plot/QC/Simulated_variants.r) to plot the calls of different methods on simulated variants.

[Codes and data](https://github.com/shishenyxx/Sperm_control_cohort_mosaicism/tree/master/Plot/Circos) for different CIRCOS plots.

-----------------------------------

### 2. Pipelines for the process of Targeted Amplicon Sequencing (TAS)
#### 2.1 Pipelines for TAS data alignment and processing

[Pipelines](https://github.com/shishenyxx/Adult_brain_somatic_mosaicism/tree/master/pipelines/MPAS_and_snMPAS_processing_pipeline) for alignment, processing, and germline variant calling of TAS reads.

#### 2.2 Pipelines for AF quantification and variant annotations

[Pipelines](https://github.com/shishenyxx/PASM/tree/master/Snakemake_pipeline) for AF quantification and variant anntations.

[Codes](https://github.com/shishenyxx/Sperm_control_cohort_mosaicism/blob/master/Pipelines/cc_validation_table_sperm_abc_saliva_blood.py) to filter and annotate on TAS data.

-----------------------------------

### 3. Pipelines for the data analysis, variant filtering, comprehensive annotations, and statistical analysis
#### 3.1 Pipelines for mosaic variant determination, annotations, and plotting

After variant calling from different strategies, variants were annotated and filtered by [a python script](https://github.com/shishenyxx/Sperm_control_cohort_mosaicism/blob/master/Pipelines/Mosaic_Calls/control_cohort_complete_from_variant_table_MSMF03.py) and positive mosaic variants as well as the corresponding samples and additional information were annotated.

[Codes](https://github.com/shishenyxx/Adult_brain_somatic_mosaicism/tree/master/permutation) for permutation analysis from gnomAD and [codes](https://github.com/shishenyxx/Sperm_control_cohort_mosaicism/tree/master/Plot/cc_annotation_counts02.py) for plotting the permutation result.

UpSet plot is generated from an [online tool](https://vcg.github.io/upset/).

#### 3.2 Pipelines for statistically analysis, and the related plotting

[Codes](https://github.com/shishenyxx/Sperm_control_cohort_mosaicism/blob/master/Mutation_accumulation_model/fit_afs.py) for the estimation of accumulation of mutations through a stepwise exponential regression regression model.

[Codes](https://github.com/shishenyxx/Sperm_control_cohort_mosaicism/blob/master/Plot/Analysis_of_validation_rates.r) for the analysis of accuracy of number of variants and estimate limit of sampling with age in different groups. 


-----------------------------------
### 4. Contact:

:email: Xiaoxu Yang: [u6055394@utah.edu](mailto:u6055394@utah.edu), [xiaoxuyanglab@gmail.com](mailto:xiaoxuyanglab@gmail.com),[xiy010@health.ucsd.edu](mailto:xiy010@health.ucsd.edu)

:email: Martin Breuss: [martin.breuss@cuanschutz.edu](mailto:martin.breuss@cuanschutz.edu)

:email: Joseph Gleeson: [jogleeson@health.ucsd.edu](mailto:jogleeson@health.ucsd.edu)

-----------------------------------
### 5. Cite the data and codes:

 <b>Yang X & Breuss MW, <i>et al.,</i> Gleeson JG. [Developmental and temporal characteristics of clonal sperm mosaicism.](https://www.sciencedirect.com/science/article/abs/pii/S0092867421008837) 2021. (<i>Cell</i>, DOI:[10.1016/j.cell.2021.07.024](https://www.doi.org/10.1016/j.cell.2021.07.024 ), PMID:[34388390](https://pubmed.ncbi.nlm.nih.gov/34388390/))</b>
 
<img src="https://user-images.githubusercontent.com/17311837/223874712-fea99e69-eb9a-4b78-a100-f4caf659b29e.png" alt="Sperm_Mosaic_Cover" width=30%> 




