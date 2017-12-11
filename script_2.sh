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


### ADD READ GROUP TO EACH ALIGNMENT
aorrg_folder='/rds/projects/2017/orsinil-bioinfoproj/VariantCall_analysis/Picard_AddOrRepleaceReadGroups/'
	mkdir -p $aorrg_folder

#### Run AddOrReplaceReadGroups to add a read group to each file. Check: https://software.broadinstitute.org/gatk/documentation/article.php?id=6472
#java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups \
#	I=$BWA_folder'Dmagna1.sorted.bam' \
#	O=$aorrg_folder'Dmagna1.ReadGroup.bam' \
#	RGID=Dmagna1 \
#	RGLB=2 \
#	RGPL=illumina \
#	RGPU=2 \
#	RGSM=Dmagna0_1
#
#java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups \
#	I=$BWA_folder'Dmagna31.sorted.bam' \
#	O=$aorrg_folder'Dmagna31.ReadGroup.bam' \
#	RGID=Dmagna31 \
#	RGLB=1 \
#	RGPL=illumina \
#	RGPU=1 \
#	RGSM=Dmagna0_1
#
#java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups \
#	I=$BWA_folder'Dmagna61.sorted.bam' \
#	O=$aorrg_folder'Dmagna61.ReadGroup.bam' \
#	RGID=Dmagna61 \
#	RGLB=2 \
#	RGPL=illumina \
#	RGPU=2 \
#	RGSM=Dmagna0_1




java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups \
	I=$BWA_folder'Dmagna12.sorted.bam' \
	O=$aorrg_folder'Dmagna12.ReadGroup.bam' \
	RGID=Dmagna12 \
	RGLB=1 \
	RGPL=illumina \
	RGPU=1 \
	RGSM=Dmagna6_2

#java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups \
#	I=$BWA_folder'Dmagna42.sorted.bam' \
#	O=$aorrg_folder'Dmagna42.ReadGroup.bam' \
#	RGID=Dmagna42 \
#	RGLB=1 \
#	RGPL=illumina \
#	RGPU=1 \
#	RGSM=Dmagna6_2
#
#java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups \
#	I=$BWA_folder'Dmagna72.sorted.bam' \
#	O=$aorrg_folder'Dmagna72.ReadGroup.bam' \
#	RGID=Dmagna72 \
#	RGLB=2 \
#	RGPL=illumina \
#	RGPU=2 \
#	RGSM=Dmagna6_2




#java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups \
#	I=$BWA_folder'Dmagna19.sorted.bam' \
#	O=$aorrg_folder'Dmagna19.ReadGroup.bam' \
#	RGID=Dmagna19 \
#	RGLB=1 \
#	RGPL=illumina \
#	RGPU=1 \
#	RGSM=Dmagna9_20
#
#java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups \
#	I=$BWA_folder'Dmagna49.sorted.bam' \
#	O=$aorrg_folder'Dmagna49.ReadGroup.bam' \
#	RGID=Dmagna49 \
#	RGLB=1 \
#	RGPL=illumina \
#	RGPU=1 \
#	RGSM=Dmagna9_20
#
#java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups \
#	I=$BWA_folder'Dmagna79.sorted.bam' \
#	O=$aorrg_folder'Dmagna79.ReadGroup.bam' \
#	RGID=Dmagna79 \
#	RGLB=2 \
#	RGPL=illumina \
#	RGPU=2 \
#	RGSM=Dmagna9_20




### MERGE WITH SAMTOOLS - Note: Samtools' merge does not reconstruct the @RG dictionary in the header. Users must provide the correct header with -h, or uses Picard which properly maintains the header dictionary in merging. 

#samtools merge -rf $merged_folder'/Dmagna-0_1.sorted.bam' $BWA_folder'Dmagna1.sorted.bam' $BWA_folder'Dmagna31.sorted.bam' $BWA_folder'Dmagna61.sorted.bam'
#samtools merge -rf $merged_folder'/Dmagna-6_2.sorted.bam' $BWA_folder'Dmagna12.sorted.bam' $BWA_folder'Dmagna42.sorted.bam' $BWA_folder'Dmagna72.sorted.bam'
#samtools merge -rf $merged_folder'/Dmagna-19_9.sorted.bam' $BWA_folder'Dmagna19.sorted.bam' $BWA_folder'Dmagna49.sorted.bam' $BWA_folder'Dmagna79.sorted.bam'


### MERGE WITH PICARD - Note that to prevent errors in downstream processing, it is critical to identify/label read groups appropriately. If different samples contain identical read group IDs, this tool will avoid collisions by modifying the read group IDs to be unique.

java -jar $EBROOTPICARD/picard.jar MergeSamFiles I=$aorrg_folder'Dmagna1.ReadGroup.bam' I=$aorrg_folder'Dmagna31.ReadGroup.bam' I=$aorrg_folder'Dmagna61.ReadGroup.bam' O=$merged_folder'/Dmagna-0_1.sorted.bam'
java -jar $EBROOTPICARD/picard.jar MergeSamFiles I=$aorrg_folder'Dmagna12.ReadGroup.bam' I=$aorrg_folder'Dmagna42.ReadGroup.bam' I=$aorrg_folder'Dmagna72.ReadGroup.bam' O=$merged_folder'/Dmagna-6_2.sorted.bam'
java -jar $EBROOTPICARD/picard.jar MergeSamFiles I=$aorrg_folder'Dmagna19.ReadGroup.bam' I=$aorrg_folder'Dmagna49.ReadGroup.bam' I=$aorrg_folder'Dmagna79.ReadGroup.bam' O=$merged_folder'/Dmagna-19_9.sorted.bam'


