#!/bin/bash

		
#SBATCH --time 8-10:0
#SBATCH --qos colboujk
#SBATCH --ntasks 10
#SBATCH --mem 50G


module purge; module load bluebear
module load apps/samtools/1.4


### Set the directories path 
home_dir='/rds/projects/2017/orsinil-bioinfoproj/CNV/'
reference=$home_dir'reference_genome/GCA_001632505.1_daphmag2.4_genomic.fna'


bam_dir=$home_dir'BAMfiles/'
merge_dir=$home_dir'merged/'

mkdir -p $merge_dir 


### Genotypes list
genotypes=$home_dir'genotypes.txt'


### Sample local name list for each genotype
samples=$home_dir'sample_mapping_list.txt'
bamList=$home_dir'bamFiles_speedseq.txt'

ls $home_dir'/BAMFiles' > bamFiles_speedseq.txt


while IFS= read -r genotype; do
	
	### Bam file merging
	
  		# Normal bam file
	samtools merge \
		-@8 \
		$merge_dir"Dmagna-"$genotype".merged.bam" \
		$bam_dir*-$genotype".bam"
	
		# Split reads bam file
  	samtools merge \
		-@8 \
		$merge_dir"Dmagna-"$genotype".merged.splitters.bam" \
		$bam_dir*-$genotype".splitters.bam"
  	
		# Discordant read pairs
	samtools merge \
		-@8 \
		$merge_dir"Dmagna-"$genotype".merged.disc.bam" \
		$bam_dir*-$genotype".discordants.bam"
  	
done < "$genotypes"




#	while IFS='	'  read localName geno; do    
#		if [ $genotype == $geno ]
#		then
#	 
#			
#		trimmed1Fq=$trim_dir$localName-$genotype-*"_1_trim.fastq.gz" 
#		trimmed2Fq=$trim_dir$localName-$genotype-*"_2_trim.fastq.gz"
#		ss_RG="@RG\tID:$localName\tSM:$genotype" 
#		ss_out=$bam_dir$localName-$genotype 
#		ss_log=$bam_dir$localName-$genotype"speedseq_align.log"
#	
#    		now=$(date +'%d-%m-%Y %r')
#    		echo " $genotype alignment is done $now"
#
#
#		### Sort and index speedseq output
#    		now=$(date +'%d-%m-%Y %r')
#    		echo "Started $genotype sorting and alignment $now"
#
#			
#			# Sort and index bam file
#		samtools sort $bam_dir$localName-$genotype".bam" \
#			-O BAM \
#			-o $sortBam_dir$localName-$genotype".sorted.bam" \
#		&& samtools index \
#			-b $sortBam_dir$localName-$genotype".sorted.bam"
#
#			### Sort and index split reads bam file
#		samtools sort $bam_dir$localName-$genotype".splitters.bam" \
#			-O BAM \
#			-o $sortBam_dir$localName-$genotype".splitters.sorted.bam" \
#		&& samtools index \
#			-b $sortBam_dir$localName-$genotype".splitters.sorted.bam"
#
#			### Sort and index discordant reads bam file 
#		samtools sort $bam_dir$localName-$genotype".discordants.bam" \
#			-O BAM \
#			-o $sortBam_dir$localName-$genotype".discordants.sorted.bam" \
#		&& samtools index \
#			-b $sortBam_dir$localName-$genotype".discordants.sorted.bam"
# 
#
#    		now=$(date +'%d-%m-%Y %r')
#    		echo " $genotype sorting and alignment is done $now"
#
#
#		fi
#
#	done < "$samples"

