#!/bin/bash

##### PART 3 
##### This script take each merged bam file and perform InDel realignment and 

### MODULE TO LOAD
   module load apps/picard/2.10.5-java-1.8.0_92    
   module load apps/gatk/v3.6

### VARIABLES FOLDERS
merged_folder='/rds/projects/2017/orsinil-bioinfoproj/VariantCall_analysis/Merged_bam_files'
project_folder='/rds/projects/2017/orsinil-bioinfoproj/VariantCall_analysis/'
Dmagna_reference=$project_folder'reference_genome/GCA_001632505.1_daphmag2.4_genomic.fna'

#### ANALYSE EACH FILE
# In order to follow all the variable substitutions, one example for each variable is added as a comment.
# All the "echo" commands are used to debug the script.

for merged_file in $project_folder*"sorted.bam"   # All the files that end with "sorted.bam" in the Merged_bam_files folder  
do

### List of the merged file
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
	sortsam_folder=$project_folder"Picard_SortSam"
	mkdir -p $sortsam_folder
	export sortsam_output=$sortsam_folder$base_name'.SortSam.bam'


		# Software commands
	java -jar $EBROOTPICARD/picard.jar SortSam I=$merged_file O=$sortsam_output
	# The option "VALIDATION STRINGENCY=SILENT" is not present in the software --help


### Run Picard MarkDuplicate on SortSam files to remove eventual PCR duplicate reads 
		
		# Initilise variables and folders 
	markduplicates_folder=$project_folder"Picard_MarkDuplicates"
	mkdir -p $markduplicates_folder
	export markduplicates_output=$markduplicates_folder$base_name'MarkDuplicates.bam'
	export markduplicates_metrics=$markduplicates_folder$base_name'MarkDuplicates.txt'

		# Software commands
	java -jar $EBROOTPICARD/picard.jar MarkDuplicates  I=$sortsam_output O=$markduplicates_output M=$markduplicates_metrics REMOVE_SEQUENCING_DUPLICATES=true 


### Run GATK RealignerTargetCreator on Markduplicate files to define intervals to target for local realignment
		
		# Initilise variables and folders 
	rtl_folder=$project_folder"GATK_RealignerTargetCreator"
	mkdir -p $rtl_folder
	export rtl_output=$rtl_folder$base_name'RealignerTargetCreator.intervals'

		# Software commands
	java -jar ${GATK_ROOT}/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $Dmagna_reference -I $markduplicates_output -o $rtl_output
	#  What is the option -L <limiting interval> ? 


### Run GATK IndelRealigner on Markduplicate files using the intevals defined by RealignerTargetCreator perform local realignment of reads around indels
		
		# Initilise variables and folders 
	indelrealigner_folder=$project_folder"GATK_IndelRealigner"
	mkdir -p $indelrealigner_folder
	export indelrealigner_output=$indelrealigner_folder$base_name'realigned.bam'

		# Software commands
	java -jar ${GATK_ROOT}/GenomeAnalysisTK.jar -T IndelRealigner -R $Dmagna_reference -I $markduplicates_output -targetIntervals $rtl_output -o $indelrealigner_output
	

done





