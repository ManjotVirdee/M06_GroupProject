#!/bin/bash

# Author:
# Gian Marco Baranzoni
		
#SBATCH --time 8-10:0
#SBATCH --qos colboujk
#SBATCH --ntasks 10
#SBATCH --mem 50G


module purge; module load bluebear
module load apps/speedseq/0.1.2
module load apps/samtools/1.4


### Set the directories path 
home_dir='/rds/projects/2017/orsinil-bioinfoproj/CNV/'
speedseq_bin=$home_dir'speedseq/bin/'
reference=$home_dir'reference_genome/GCA_001632505.1_daphmag2.4_genomic.fna'


trim_dir=$home_dir'trimReads/'
bam_dir=$home_dir'bam/'
sortBam_dir=$home_dir'sortBam/'
merge_dir=$home_dir'merged/'
temp_dir=$home_dir'temp'

mkdir -p $bam_dir $sortBam_dir $merge_dir $temp_dir 


### Genotypes list
genotypes=$home_dir'genotypes.txt'


### Sample local name list for each genotype
samples=$home_dir'sample_mapping_list.txt'


while IFS= read -r genotype; do

	while IFS='	'  read localName geno; do    
		
#			echo "$genotype"
#			echo "$localName"
#			echo "$geno"
#			echo ""
		
		if [ $genotype == $geno ]
		then
	 
			
		### Speedseq align
    		now=$(date +'%d-%m-%Y %r')
    		echo "Started $genotype align $now" 
    		
		trimmed1Fq=$trim_dir$localName-$genotype-*"_1_trim.fastq.gz" 
		trimmed2Fq=$trim_dir$localName-$genotype-*"_2_trim.fastq.gz"
		ss_RG="@RG\tID:$localName\tSM:$genotype" 
		ss_out=$bam_dir$localName-$genotype 
		ss_log=$bam_dir$localName-$genotype"speedseq_align.log"
	
#			echo $trimmed1Fq
#			echo $trimmed2Fq
#			echo $ss_RG
#			echo $ss_out
#			echo $ss_log
	
#			speedseq align\
		bash /rds/projects/2017/orsinil-bioinfoproj/CNV/speedseq/bin/speedseq align \
			-t 10 \
			-T $temp_dir \
			-o $ss_out \
			-R $ss_RG \
			$reference \
			$trimmed1Fq \
			$trimmed2Fq 

#&> $ss_log

    		now=$(date +'%d-%m-%Y %r')
    		echo " $genotype alignment is done $now"


		### Sort and index speedseq output
    		now=$(date +'%d-%m-%Y %r')
    		echo "Started $genotype sorting and alignment $now"

			
			# Sort and index bam file
		samtools sort $bam_dir$localName-$genotype".bam" \
			-O BAM \
			-o $sortBam_dir$localName-$genotype".sorted.bam" \
		&& samtools index \
			-b $sortBam_dir$localName-$genotype".sorted.bam"

			### Sort and index split reads bam file
		samtools sort $bam_dir$localName-$genotype".splitters.bam" \
			-O BAM \
			-o $sortBam_dir$localName-$genotype".splitters.sorted.bam" \
		&& samtools index \
			-b $sortBam_dir$localName-$genotype".splitters.sorted.bam"

			### Sort and index discordant reads bam file 
		samtools sort $bam_dir$localName-$genotype".discordants.bam" \
			-O BAM \
			-o $sortBam_dir$localName-$genotype".discordants.sorted.bam" \
		&& samtools index \
			-b $sortBam_dir$localName-$genotype".discordants.sorted.bam"
 

    		now=$(date +'%d-%m-%Y %r')
    		echo " $genotype sorting and alignment is done $now"


		fi

	done < "$samples"

   
	now=$(date +'%d-%m-%Y %r')
	echo " Alignment for all individuals is done $now"

	
	### Bam file merging
	now=$(date +'%d-%m-%Y %r')
   	echo "Merging all the bam files is done $now"

	
  		# Normal bam file
	samtools merge \
		-@8 \
		$merge_dir$genotype".merged.bam" \
		$sortBam_dir*-$genotype".bam"
	
		# Split reads bam file
  	samtools merge \
		-@8 \
		$merge_dir$genotype".merged.splitters.bam" \
		$sortBam_dir*-$genotype".splitters.bam"
  	
		# Discordant read pairs
	samtools merge \
		-@8 \
		$merge_dir$genotype".merged.disc.bam" \
		$sortBam_dir*-$genotype".discordants.bam"
  	
	now=$(date +'%d-%m-%Y %r')
   	echo "Merging all the bam files is done $now"

done < "$genotypes"





