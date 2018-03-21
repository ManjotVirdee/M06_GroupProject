#!bin/bash

Authors:
Gian Marco Baranzoni and Manjot Virdee
 
project_folder='/rds/projects/2017/orsinil-bioinfoproj/'
vcf_file_orig=$project_folder'Dmagna_30_individuals_SNPs.vcf.gz'
gwas_folder='/rds/homes/g/gxb738/gwas_plink/'
vcf_file=$gwas_folder'Dmagna_30_individuals_SNPs.vcf'
vcfCD_file=$project_folder'vcfAnn/Dmagna30_CD.vcf'





#### Run Plink
##module load apps/plink/1.07
#
## --double-id  -  both family and within-family IDs to be set to the sample ID
## --allow-extra-chr  -  to use scaffolds as cromosome ID 
## --geno - filters out all variants with missing call rates exceeding the provided value (default 0.1) to be removed. I choose 0.03 to excluded all the snps with a missing value in any samples
## --maf filters out all variants with minor allele frequency below the provided threshold (default 0.01), while --max-maf imposes an upper MAF bound. Reed et al (2015) suggest 5% for small samples, which corresponds to 1.65 individuals. Since our number of samples is very low, I will start from 10% (3.3 individuals) and we can reduce it later.
## --hwe filters out all variants which have Hardy-Weinberg equilibrium exact test p-value below the provided threshold. We recommend setting a low threshold—serious genotyping errors often yield extreme p-values like 1e-50 which are detected by any reasonable configuration of this test, while genuine SNP-trait associations can be expected to deviate slightly from Hardy-Weinberg equilibrium (so it's dangerous to choose a threshold that filters out too many variants).


## Old run with wrong file 
#plink\
#	--vcf $vcf_file \
#	--double-id \
#	--allow-extra-chr \
#	--allow-no-sex \
#	--geno 0.03 \
#	--maf 0.1 \
#	--hwe 1e-50 \
#	--pheno temp_37.pheno \
#	--all-pheno \
#	--assoc \
#	--linear \
#	--out temp_37_out



##### Run plink with phenotype data for each temperature
#
#plink\
#	--vcf $vcf_file \
#	--double-id \
#	--allow-extra-chr \
#	--allow-no-sex \
#	--geno 0.03 \
#	--maf 0.1 \
#	--hwe 1e-50 \
#	--pheno temp_18.pheno \
#	--all-pheno \
#	--assoc \
#	--out temp_18/temp_18_out
#
#plink\
#	--vcf $vcf_file \
#	--double-id \
#	--allow-extra-chr \
#	--allow-no-sex \
#	--geno 0.03 \
#	--maf 0.1 \
#	--hwe 1e-50 \
#	--pheno temp_21.pheno \
#	--all-pheno \
#	--assoc \
#	--out temp_21/temp_21_out
#
#plink\
#	--vcf $vcf_file \
#	--double-id \
#	--allow-extra-chr \
#	--allow-no-sex \
#	--geno 0.03 \
#	--maf 0.1 \
#	--hwe 1e-50 \
#	--pheno temp_24.pheno \
#	--all-pheno \
#	--assoc \
#	--out temp_24/temp_24_out


#
##### Run plink with phenotype data with the difference between the temparatures. Columns of population and mortality were removed.
#
#plink\
#	--vcf $vcf_file \
#	--double-id \
#	--allow-extra-chr \
#	--allow-no-sex \
#	--geno 0.03 \
#	--maf 0.1 \
#	--hwe 1e-50 \
#	--pheno temp_T18d21.pheno \
#	--all-pheno \
#	--assoc \
#	--out temp_T18d21/temp_T18d21_out
#
#plink\
#	--vcf $vcf_file \
#	--double-id \
#	--allow-extra-chr \
#	--allow-no-sex \
#	--geno 0.03 \
#	--maf 0.1 \
#	--hwe 1e-50 \
#	--pheno temp_T21d24.pheno \
#	--all-pheno \
#	--assoc \
#	--out temp_T21d24/temp_T21d24_out
#
#plink\
#	--vcf $vcf_file \
#	--double-id \
#	--allow-extra-chr \
#	--allow-no-sex \
#	--geno 0.03 \
#	--maf 0.1 \
#	--hwe 1e-50 \
#	--pheno temp_T18d24.pheno \
#	--all-pheno \
#	--assoc \
#	--out temp_T18d24/temp_T18d24_out
#
#
#
#
#
#
#
##### Run plink with vcf file with variants only in the coding regions

plink\
	--vcf $vcfCD_file \
	--double-id \
	--allow-extra-chr \
	--allow-no-sex \
	--geno 0.03 \
	--maf 0.1 \
	--hwe 1e-50 \
	--pheno temp_18.pheno \
	--all-pheno \
	--assoc \
	--out temp_18CD/temp_18CD_out


plink\
	--vcf $vcfCD_file \
	--double-id \
	--allow-extra-chr \
	--allow-no-sex \
	--geno 0.03 \
	--maf 0.1 \
	--hwe 1e-50 \
	--pheno temp_T18d24.pheno \
	--all-pheno \
	--assoc \
	--out temp_T18d24CD/temp_T18d24CD_out

