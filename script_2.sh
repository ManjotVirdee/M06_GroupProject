#!/bin/bash

##### PART 2 
##### This script take the three technical replicates in a sorted bam format and merge them in a single file.

### MODULE TO LOAD
#   module load apps/samtools
#   module load apps/picard/2.10.5-java-1.8.0_92


### VARIABLES FOLDERS
BWA_folder='/rds/projects/2017/orsinil-bioinfoproj/VariantCall_analysis/BWA/'
merged_folder='/rds/projects/2017/orsinil-bioinfoproj/VariantCall_analysis/Merged_bam_files'
	mkdir -p $merged_folder

### MERGE WITH SAMTOOLS - Note: Samtools' merge does not reconstruct the @RG dictionary in the header. Users must provide the correct header with -h, or uses Picard which properly maintains the header dictionary in merging. 

#samtools merge -rf $merged_folder'/Dmagna-0_1.sorted.bam' $BWA_folder'Dmagna1.sorted.bam' $BWA_folder'Dmagna31.sorted.bam' $BWA_folder'Dmagna61.sorted.bam'
#samtools merge -rf $merged_folder'/Dmagna-6_2.sorted.bam' $BWA_folder'Dmagna12.sorted.bam' $BWA_folder'Dmagna42.sorted.bam' $BWA_folder'Dmagna72.sorted.bam'
#samtools merge -rf $merged_folder'/Dmagna-19_9.sorted.bam' $BWA_folder'Dmagna19.sorted.bam' $BWA_folder'Dmagna49.sorted.bam' $BWA_folder'Dmagna79.sorted.bam'


### MERGE WITH PICARD - Note that to prevent errors in downstream processing, it is critical to identify/label read groups appropriately. If different samples contain identical read group IDs, this tool will avoid collisions by modifying the read group IDs to be unique.

java -jar $EBROOTPICARD/picard.jar MergeSamFiles I=$BWA_folder'Dmagna1.sorted.bam' I=$BWA_folder'Dmagna31.sorted.bam' I=$BWA_folder'Dmagna61.sorted.bam' O=$merged_folder'/Dmagna-0_1.sorted.bam'
java -jar $EBROOTPICARD/picard.jar MergeSamFiles I=$BWA_folder'Dmagna12.sorted.bam' I=$BWA_folder'Dmagna42.sorted.bam' I=$BWA_folder'Dmagna72.sorted.bam' O=$merged_folder'/Dmagna-6_2.sorted.bam'
java -jar $EBROOTPICARD/picard.jar MergeSamFiles I=$BWA_folder'Dmagna19.sorted.bam' I=$BWA_folder'Dmagna49.sorted.bam' I=$BWA_folder'Dmagna79.sorted.bam' O=$merged_folder'/Dmagna-19_9.sorted.bam'


