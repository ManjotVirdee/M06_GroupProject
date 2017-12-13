#!/bin/bash

##### PART 1 
##### Quality control and aligment of Daphina magna genomes in the "DmagnaSamples_subset" directory


### MODULE TO LOAD
#   module load apps/fastqc/v0.11.2  
#   module load apps/trimmomatic
#   module load apps/samtools
#   module load apps/BWA


### VARIABLES FOLDERS
project_folder='/rds/projects/2017/orsinil-bioinfoproj/VariantCall_analysis/'
raw_data_folder=$project_folder'DmagnaSamples_subset/'


### Indexing the reference genome (bwtsw for large genomes)

	# Initilise variables 
export Dmagna_reference=$project_folder'reference_genome/GCA_001632505.1_daphmag2.4_genomic.fna'

	# Software commands
# bwa index -a bwtsw $Dmagna_reference  # bwtsw for large genomes



#### ANALYSE EACH FILE
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
	export fastq_base="${fastq_file_forward%%-*}"
		# Dmagna1
		echo $fastq_forward
		echo $fastq_reverse
		echo $fastq_file_forward
		echo $fastq_file_reverse
		echo $fastq_base


### Run FastQC on raw data
		
		# Initilise variables and folders 
	fastqc_folder=$project_folder"FastQC"
	mkdir -p $fastqc_folder
	fastqc_rawdata=$project_folder"FastQC/raw_data"
	mkdir -p $fastqc_rawdata
		
		# Software commands
	fastqc $fastq_forward -o $fastqc_rawdata
	fastqc $fastq_reverse -o $fastqc_rawdata


### Run Trimmomatic on raw data

		# Initilise variables and folders 
	trimmomatic_data=$project_folder"Trimmomatic_data/"
	mkdir -p $trimmomatic_data
	export trimmed_output=$trimmomatic_data$fastq_base
	export adapter_fasta=$project_folder"TruSeq3-PE-2.fa"

		# Software commands
	trimmomatic PE $fastq_forward $fastq_reverse -baseout $trimmed_output CROP:225 ILLUMINACLIP:$adapter_fasta:2:30:10 LEADING:10 TRAILING:10 SLIDINGWINDOW:50:20 MINLEN:36
	

### Run FastQC of Trimmomatic data

		# Initilise variables and folders 
	fastqc_trimmomatic=$project_folder"FastQC/Trimmomatic_data"
	mkdir -p $fastqc_trimmomatic

		# Software commands
	fastqc $trimmed_output"_1P" -o $fastqc_trimmomatic
	fastqc $trimmed_output"_2P" -o $fastqc_trimmomatic
	fastqc $trimmed_output"_1U" -o $fastqc_trimmomatic
	fastqc $trimmed_output"_2U" -o $fastqc_trimmomatic


### Run BWA of Trimmomatic data

		# Initilise variables and folders 
	bwa_data=$project_folder"BWA/"
	mkdir -p $bwa_data
	export bwa_output=$bwa_data$fastq_base".sorted.bam"

		# Software commands
	bwa mem $Dmagna_reference $trimmed_output"_1P" $trimmed_output"_2P" | samtools view -bh | samtools sort > $bwa_output

	# I should have added the option -M to Mark shorter split hits as secondary (for Picard compatibility)

done





