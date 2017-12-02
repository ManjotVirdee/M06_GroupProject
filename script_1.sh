#!/bin/bash

##### FILTERING AND QUALITY CONTROL
### MODULE TO LOAD
#   module load apps/fastqc/v0.11.2  
#   module load apps/trimgalore/v0.4.2
#   module load apps/trimmomatic


### VARIABLES FOLDERS
project_folder='/rds/projects/2017/orsinil-bioinfoproj/VariantCall_analysis/'
raw_data_folder=$project_folder'DmagnaSamples_subset/'

### MAKE DIRECTORIES - the option p is to avoid the creation of a directory if it already exist
fastqc_folder=$project_folder"FastQC"
	mkdir -p $fastqc_folder
fastqc_rawdata=$project_folder"FastQC/raw_data"
	mkdir -p $fastqc_raw_data
fastqc_trimgalore=$project_folder"FastQC/TrimGalore_data"
	mkdir -p $fastqc_trimgalore
fastqc_trimmomatic=$project_folder"FastQC/Trimmomatic_data"
	mkdir -p $project_trimmomatic
trimgalore_data=$project_folder"TrimGalore_data/"
	mkdir -p $trimgalore_data
trimmomatic_data=$project_folder"Trimmomatic_data/"
	mkdir -p $trimmomatic_data


#### RUN EACH FILE
# In order to follow all the variable substitutions, one example for each variable is added as a comment.
# All the "echo" commands are used to debug the script.

for fastq_forward in $raw_data_folder*"-R1.fastq.gz";   # All the files that end with "-R1.fastq.gz"  
do
	### List of the Daphnias
	export fastq_forward 
		# /rds/projects/2017/orsinil-bioinfoproj/VariantCall_analysis/DmagnaSamples_subset/Dmagna1-0_1-L002-R1.fastq.gz
	export fastq_reverse="$( echo $fastq_forward | sed -e 's/-R1.fastq.gz/-R2.fastq.gz/')" 
		# /rds/projects/2017/orsinil-bioinfoproj/VariantCall_analysis/DmagnaSamples_subset/Dmagna1-0_1-L002-R2.fastq.gz
	export fastq_file_forward="$( echo $fastq_forward | sed -e 's/\/rds\/projects\/2017\/orsinil-bioinfoproj\/VariantCall_analysis\/DmagnaSamples_subset\///')" 
		# Dmagna1-0_1-L002-R1.fastq.gz
	export fastq_file_reverse="$( echo $fastq_file_forward | sed -e 's/-R1.fastq.gz/-R2.fastq.gz/')" 
		# Dmagna1-0_1-L002-R2.fastq.gz
	export fastq_base="${fastq_file_forward%%-*}" # I need to change this!
		# Dmagna1 
		echo $fastq_forward
		echo $fastq_reverse
		echo $fastq_file_forward
		echo $fastq_file_reverse
		echo $fastq_base


	# Run FastQC on raw data
	fastqc $fastq_forward -o $fastqc_rawdata
	fastqc $fastq_reverse -o $fastqc_rawdata

	# Run TrimGalore on raw data
	trim_galore --trim-n --illumina --paired $fastq_forward $fastq_reverse -o $trimgalore_data
	mv $trimgalore_data$fastq_file_forward"_trimming_report.txt" $fastqc_trimgalore
	mv $trimgalore_data$fastq_file_reverse"_trimming_report.txt" $fastqc_trimgalore

	# Run FastQC on TrimGalore data
	fastqc $fastq_forward -o $fastqc_trimgalore
	fastqc $fastq_reverse -o $fastqc_trimgalore

	# Run Trimmomatic on raw data
	export trimgalore_forward=$trimgalore_data"$( echo $fastq_file_forward | sed -e 's/.fastq.gz/_val_1.fq.gz/')"
	export trimgalore_reverse=$trimgalore_data"$( echo $fastq_file_reverse | sed -e 's/.fastq.gz/_val_2.fq.gz/')"
	export trimmed_output=$trimmomatic_data$fastq_base
	trimmomatic PE -phred33 $trimgalore_forward $trimgalore_reverse -baseout $trimmed_output LEADING:3 TRAILING:3 SLIDINGWINDOW:80:25 MINLEN:35 

	echo $trimgalore_forward
	echo $trimgalore_reverse
	echo $trimmed_output

	# Run FastQC of Trimmomatic data
	fastqc $trimmed_output"_1P" -o $fastqc_trimmomatic
	fastqc $trimmed_output"_2P" -o $fastqc_trimmomatic

echo ""
done





