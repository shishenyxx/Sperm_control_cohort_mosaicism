##This scripts contains the workflow for whole-genome variant extraction and PCA, codes from Danny Antaki, modified by An Nguyen and Xiaoxu Yang
#First generate gvcf from haplotypecaller, use -D <dbSNP_file>
cd Ctrl_cohort_gvcf/gvcfs
1. Combine gvcfs
vcf-merge ID01-Blood.vcf.gz ID02-Blood.vcf.gz ID03-Blood.vcf.gz ID04-Blood.vcf.gz ID05-Blood.vcf.gz ID06-Blood.vcf.gz ID07-Blood.vcf.gz ID08-Blood.vcf.gz ID09-Blood.vcf.gz ID10-Blood.vcf.gz ID11-Blood.vcf.gz ID12-Blood.vcf.gz ID13-Blood.vcf.gz ID14-Blood.vcf.gz ID15-Blood.vcf.gz ID16-Blood.vcf.gz ID17-Blood.vcf.gz|bgzip -c >Control_cohort_blood.vcf.gz

2. reformat your VCF
bcftools view -i'TYPE=="SNP" && N_ALT=1' Control_cohort_blood.vcf.gz | bcftools annotate -Oz -x ID -I +'%CHROM:%POS:%REF:%ALT' >in.snv.reid.genotypes.vcf.gz

3. convert VCF to plink bfiles
plink --vcf in.snv.reid.genotypes.vcf.gz --make-bed --out my_bfile --double-id

4. Extract the markers (remove markers in the repetitive regions)
my_join.pl -F 2 -f 2 -a my_bfile.bim -b 1kgp.pruned.bim -m|awk '{print $2}'|sort|uniq > plink.prune.in

cat plink.prune.in|sed 's/:/\t/g'|awk -v OFS="\t" '{print $1,$2,$2,$3,$4}'>plink.prune.bed

bedtools annotate -counts -i plink.prune.bed -files resources/hg19_lumpy_exclude_merged.bed.gz|awk '{if ($6=="0") print $1":"$2":"$4":"$5}'>plink.prune.in

5. extract pruned markers from combined vcf
plink --bfile my_bfile --extract plink.prune.in --make-bed --out my_bfile.pruned

6. merge 1kgp with combined vcf
plink --bfile 1kgp.pruned --bmerge my_bfile.pruned --geno 0 --merge-mode 1 --make-bed --out 1kgp_control_cohort

7. pca (10 components)
plink --bfile 1kgp_control_cohort --pca 10 --out control_cohort.pca

