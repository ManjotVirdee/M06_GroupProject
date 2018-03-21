### VariantAnnotation - Server ###

# Author:
# Gian Marco Baranzoni

### Libraries
library(VariantAnnotation)
library(GenomicFeatures)
library(Rsamtools)

setwd("/rds/projects/2017/orsinil-bioinfoproj/vcfAnn")

## import a vcf file
# Index the vcf using tabix in unix using SAMtools as follow: "tabix -p vcf myfile.vcf"
DMvcf_tab <- TabixFile("/rds/projects/2017/orsinil-bioinfoproj/vcfAnn/Dmagna_30_individuals_SNPs.vcf.gz")

DMvcf <-  readVcf ("/rds/projects/2017/orsinil-bioinfoproj/vcfAnn/Dmagna_30_individuals_SNPs.vcf.gz") 

## import the annotation as GRange object
DMtxdb_scaffID <- makeTxDbFromGFF(file="/rds/projects/2017/orsinil-bioinfoproj/vcfAnn/Dmagna_scaffID.gff", format="auto")

## extract variants in coding regions
GRanges_cr <- locateVariants(DMvcf, DMtxdb_scaffID, CodingVariants())

DMvcf_cr <-  readVcf ("/rds/projects/2017/orsinil-bioinfoproj/vcfAnn/Dmagna_30_individuals_SNPs.vcf.gz", seqinfo(DMvcf), GRanges_cr) 

writeVcf(DMvcf_cr, "/rds/projects/2017/orsinil-bioinfoproj/vcfAnn/Dmagna30_CD.vcf")
