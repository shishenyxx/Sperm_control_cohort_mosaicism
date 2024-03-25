# use this script to help format the output, it will identify population groups in 1kgp phase 3

/home/dantakli/4dbsm/ancestry/src/tabulate_pca.pl

1. reformat your VCF

bcftools view -i'TYPE=="SNP" && N_ALT=1' in.vcf.gz | bcftools annotate -Oz -x ID -I +'%CHROM:%POS:%REF:%ALT' >in.snv.reid.genotypes.vcf.gz

2. convert VCF to plink bfiles

plink --vcf in.snv.reid.vcf.gz --make-bed --out my_bfile --double-id

3. extract pruned markers from your vcf

plink --bfile my_bfile --extract /projects/ps-gleesonlab6/4dbsm/jun_2019/hpcaller/ancestry/plink.prune.in --make-bed --out my_bfile.pruned

5. merge 1kgp with your vcf

plink --bfile 1kgp.pruned --bmerge my_bfile.pruned --geno 0 --merge-mode 1 --make-bed --out 1kgp_my_bfile

6. pca (10 components)

plink --bfile 1kgp_my_bfile --pca --out my_bfile.pca

(plinkv1.9 from Danny: /home/dantakli/4dbsm/ancestry/README.md)

(Also refer to /projects/ps-gleesonlab6/4dbsm/jun_2019/hpcaller/ancestry/)

# these are the LD pruned markers

/projects/ps-gleesonlab6/4dbsm/jun_2019/hpcaller/ancestry/1kgp.pruned.*

#PCA command

plink --bfile $bfile_prefix --pca --out $bfile\.pca
