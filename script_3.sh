#!/bin/bash

##### PART 3 
##### This script take each merged bam file and perform InDel realignment and 

### MODULE TO LOAD
#   module load apps/picard/2.10.5-java-1.8.0_92    
#   module load apps/gatk/v3.6

### VARIABLES FOLDERS
merged_folder='/rds/projects/2017/orsinil-bioinfoproj/VariantCall_analysis/Merged_bam_files/'
project_folder='/rds/projects/2017/orsinil-bioinfoproj/VariantCall_analysis/'
Dmagna_reference=$project_folder'reference_genome/GCA_001632505.1_daphmag2.4_genomic.fna'


### CREATE REFERENCE DICTIONARY
# java -jar $EBROOTPICARD/picard.jar CreateSequenceDictionary R=$Dmagna_reference O=$Dmagna_reference'.dict'


#### ANALYSE EACH FILE
# In order to follow all the variable substitutions, one example for each variable is added as a comment.
# All the "echo" commands are used to debug the script.

for merged_file in $merged_folder*"sorted.bam"   # All the files that end with "sorted.bam" in the Merged_bam_files folder  
do

#### List of the merged file
	export merged_file 
		#  /rds/projects/2017/orsinil-bioinfoproj/VariantCall_analysis/Merged_bam_files/Dmagna-0_1.sorted.bam
	export temp_name="$( echo $merged_file | sed -e 's/\/rds\/projects\/2017\/orsinil-bioinfoproj\/VariantCall_analysis\/Merged_bam_files\///')" 
		# Dmagna-0_1.sorted.bam
	export base_name="${temp_name%%.*}"
		# Dmagna-0_1
		echo $merged_file
		echo $temp_name
		echo $base_name


### Run Picard SortSam on merged files to re-sort the reads 
		
		# Initilise variables and folders 
	ss_folder=$project_folder"Picard_SortSam/"
	mkdir -p $ss_folder
	export ss_output=$ss_folder$base_name'.SortSam.bam'

		# Software commands
	# not working
	# java -jar $EBROOTPICARD/picard.jar SortSam I=$merged_file O=$ss_output VALIDATION STRINGENCY=SILENT 
	# corrected
#	java -jar $EBROOTPICARD/picard.jar SortSam I=$merged_file O=$ss_output VALIDATION_STRINGENCY=SILENT SORT_ORDER=coordinate CREATE_INDEX=true

	
### Run GATK DepthOfCoverage to asses sequence coverage of SortSam files  
		
		# Initilise variables and folders 
	doc_folder=$project_folder"DepthOfCoverage/"
	doc_ss_folder=$doc_folder"SortSam_files/"
	mkdir -p $doc_folder
	mkdir -p $doc_ss_folder
	export doc_ss_output=$doc_ss_folder$base_name

		echo ""	
		echo $doc_folder
		echo $doc_ss_folder
		echo $doc_ss_output
		echo $Dmagna_reference
		
		# Software commands
#	java -jar ${GATK_ROOT}/GenomeAnalysisTK.jar -T DepthOfCoverage -R $Dmagna_reference -I $ss_output -o $doc_ss_output
	

### Run Picard MarkDuplicate on SortSam files to remove eventual PCR duplicate reads 
		
		# Initilise variables and folders 
	md_folder=$project_folder"Picard_MarkDuplicates/"
	mkdir -p $md_folder
	#mkdir -p $md_folder"temp"
	export md_output=$md_folder$base_name'MarkDuplicates.bam'
	export md_metrics=$md_folder$base_name'MarkDuplicates.txt'

		echo $ss_output

		# Software commands 
	java -Xmx8G -jar $EBROOTPICARD/picard.jar MarkDuplicates I=$ss_output O=$md_output M=$md_metrics ASSUME_SORT_ORDER=coordinate REMOVE_SEQUENCING_DUPLICATES=true VALIDATION_STRINGENCY=SILENT TMP_DIR=/temp


### Run GATK RealignerTargetCreator on Markduplicate files to define intervals to target for local realignment
		
		# Initilise variables and folders 
	rtl_folder=$project_folder"GATK_RealignerTargetCreator/"
	mkdir -p $rtl_folder
	export rtl_output=$rtl_folder$base_name'RealignerTargetCreator.intervals'

		# Software commands
#	java -jar ${GATK_ROOT}/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $Dmagna_reference -I $ss_output -o $rtl_output


### Run GATK IndelRealigner on Markduplicate files using the intevals defined by RealignerTargetCreator perform local realignment of reads around indels
		
		# Initilise variables and folders 
	ir_folder=$project_folder"GATK_IndelRealigner/"
	mkdir -p $ir_folder
	export ir_output=$ir_folder$base_name'realigned.bam'

		# Software commands
#	java -jar ${GATK_ROOT}/GenomeAnalysisTK.jar -T IndelRealigner -R $Dmagna_reference -I $ss_output -targetIntervals $rtl_output -o $ir_output
	

### Run GATK DepthOfCoverage to asses sequence coverage of IndelRealigner files  
		
		# Initilise variables and folders 
	doc_ir_folder=$doc_folder"IndelRealigner_files/"
	mkdir -p $doc_ir_folder
	export doc_ir_output=$doc_ir_folder$base_name

		# Software commands
#	java -jar ${GATK_ROOT}/GenomeAnalysisTK.jar -T DepthOfCoverage -R $Dmagna_reference -I $ir_output -o $doc_ir_output


### Run Samtools mpileup on GATK IndelRealigner files
		
		# Initilise variables and folders 
	mpileup_folder=$project_folder"mpileup/"
	mkdir -p $mpileup_folder
	export mpileup_output=$mpileup_folder$base_name'.vcf.gz'

		# Software commands
#	samtools mpileup -f $Dmagna_reference $ir_output -v > $mpileup_output


done




