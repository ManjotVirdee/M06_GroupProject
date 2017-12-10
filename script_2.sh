#!/bin/bash

##### PART 2 
##### This script take the three technical replicates in a sorted bam format and merge them in a single file.

### MODULE TO LOAD
#   module load apps/samtools


### VARIABLES FOLDERS
BWA_folder='/rds/projects/2017/orsinil-bioinfoproj/VariantCall_analysis/BWA/'
merged_folder='/rds/projects/2017/orsinil-bioinfoproj/VariantCall_analysis/Merged_bam_files'
	mkdir -p $merged_folder

### MERGE WITH SAMTOOLS 
samtools merge -f $merged_folder'/Dmagna-0_1.sorted.bam' $BWA_folder'Dmagna1.sorted.bam' $BWA_folder'Dmagna31.sorted.bam' $BWA_folder'Dmagna61.sorted.bam'
samtools merge -f $merged_folder'/Dmagna-6_2.sorted.bam' $BWA_folder'Dmagna12.sorted.bam' $BWA_folder'Dmagna42.sorted.bam' $BWA_folder'Dmagna72.sorted.bam'
samtools merge -f $merged_folder'/Dmagna-19_9.sorted.bam' $BWA_folder'Dmagna19.sorted.bam' $BWA_folder'Dmagna49.sorted.bam' $BWA_folder'Dmagna79.sorted.bam'

