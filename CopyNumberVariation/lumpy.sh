#!/bin/bash

		
#SBATCH --time 8-10:0
#SBATCH --qos colboujk
#SBATCH --ntasks 5
#SBATCH --mem 50G


#module purge; module load bluebear
#module load apps/lumpy/v0.2.12

### list of bam files were generated as follow:
#ls BAMFiles/ | grep sorted.merged.bam | sort >> lumpy.sh 
#ls BAMFiles/ | grep splitters.merged.bam | sort >> lumpy.sh 
#ls BAMFiles/ | grep discordants.sorted.merged.bam | sort >> lumpy.sh 

### Set the directories path 
home_dir='/rds/projects/2017/orsinil-bioinfoproj/CNV/'
temp_dir=$home_dir'temp/'
lumExp_dir=$home_dir'LumpyExpress/'
bam_dir=$home_dir'BAMFiles/'

mkdir -p $lumExp_dir

### Run Lumpy-express
/gpfs/bb/dhandapv/local/miniconda2/bin/lumpyexpress \
	-B $bam_dir"0_1.sorted.merged.bam",\
		$bam_dir"0_2.sorted.merged.bam",\
		$bam_dir"0_4.sorted.merged.bam",\
		$bam_dir"12_3.sorted.merged.bam",\
		$bam_dir"12_4.sorted.merged.bam",\
		$bam_dir"12_5_1.sorted.merged.bam",\
		$bam_dir"1_2.sorted.merged.bam",\
		$bam_dir"13_1.sorted.merged.bam",\
		$bam_dir"13_2.sorted.merged.bam",\
		$bam_dir"13_3.sorted.merged.bam",\
		$bam_dir"13_5_1.sorted.merged.bam",\
		$bam_dir"14_5_1.sorted.merged.bam",\
		$bam_dir"15_5_1.sorted.merged.bam",\
		$bam_dir"16_5_1.sorted.merged.bam",\
		$bam_dir"16_5_2.sorted.merged.bam",\
		$bam_dir"17_5_1.sorted.merged.bam",\
		$bam_dir"17_5_2.sorted.merged.bam",\
		$bam_dir"18_5A.sorted.merged.bam",\
		$bam_dir"19_1B.sorted.merged.bam",\
		$bam_dir"19_3.sorted.merged.bam",\
		$bam_dir"19_4.sorted.merged.bam",\
		$bam_dir"19_5_3.sorted.merged.bam",\
		$bam_dir"19_5.sorted.merged.bam",\
		$bam_dir"20_5.sorted.merged.bam",\
		$bam_dir"21_5_1.sorted.merged.bam",\
		$bam_dir"21_5_2.sorted.merged.bam",\
		$bam_dir"2_1.sorted.merged.bam",\
		$bam_dir"22_5_2B.sorted.merged.bam",\
		$bam_dir"2_5_11.sorted.merged.bam",\
		$bam_dir"2_5_9.sorted.merged.bam",\
		$bam_dir"3_5_15.sorted.merged.bam",\
		$bam_dir"3_5_1.sorted.merged.bam",\
		$bam_dir"3_5_2.sorted.merged.bam",\
		$bam_dir"3_6.sorted.merged.bam",\
		$bam_dir"6_2.sorted.merged.bam",\
		$bam_dir"6_3.sorted.merged.bam",\
		$bam_dir"7_3.sorted.merged.bam",\
		$bam_dir"7_5_4.sorted.merged.bam",\
		$bam_dir"7_5.sorted.merged.bam",\
		$bam_dir"8_5_3.sorted.merged.bam",\
		$bam_dir"9_20.sorted.merged.bam",\
		$bam_dir"9_5_1.sorted.merged.bam",\
		$bam_dir"9_5_3.sorted.merged.bam",\
		$bam_dir"9_6.sorted.merged.bam",\
	-S $bam_dir"0_1.splitters.merged.bam",\
		$bam_dir"0_2.splitters.merged.bam",\
		$bam_dir"0_4.splitters.merged.bam",\
		$bam_dir"12_3.splitters.merged.bam",\
		$bam_dir"12_4.splitters.merged.bam",\
		$bam_dir"12_5_1.splitters.merged.bam",\
		$bam_dir"1_2.splitters.merged.bam",\
		$bam_dir"13_1.splitters.merged.bam",\
		$bam_dir"13_2.splitters.merged.bam",\
		$bam_dir"13_3.splitters.merged.bam",\
		$bam_dir"13_5_1.splitters.merged.bam",\
		$bam_dir"14_5_1.splitters.merged.bam",\
		$bam_dir"15_5_1.splitters.merged.bam",\
		$bam_dir"16_5_1.splitters.merged.bam",\
		$bam_dir"16_5_2.splitters.merged.bam",\
		$bam_dir"17_5_1.splitters.merged.bam",\
		$bam_dir"17_5_2.splitters.merged.bam",\
		$bam_dir"18_5A.splitters.merged.bam",\
		$bam_dir"19_1B.splitters.merged.bam",\
		$bam_dir"19_3.splitters.merged.bam",\
		$bam_dir"19_4.splitters.merged.bam",\
		$bam_dir"19_5_3.splitters.merged.bam",\
		$bam_dir"19_5.splitters.merged.bam",\
		$bam_dir"20_5.splitters.merged.bam",\
		$bam_dir"21_5_1.splitters.merged.bam",\
		$bam_dir"21_5_2.splitters.merged.bam",\
		$bam_dir"2_1.splitters.merged.bam",\
		$bam_dir"22_5_2B.splitters.merged.bam",\
		$bam_dir"2_5_11.splitters.merged.bam",\
		$bam_dir"2_5_9.splitters.merged.bam",\
		$bam_dir"3_5_15.splitters.merged.bam",\
		$bam_dir"3_5_1.splitters.merged.bam",\
		$bam_dir"3_5_2.splitters.merged.bam",\
		$bam_dir"3_6.splitters.merged.bam",\
		$bam_dir"6_2.splitters.merged.bam",\
		$bam_dir"6_3.splitters.merged.bam",\
		$bam_dir"7_3.splitters.merged.bam",\
		$bam_dir"7_5_4.splitters.merged.bam",\
		$bam_dir"7_5.splitters.merged.bam",\
		$bam_dir"8_5_3.splitters.merged.bam",\
		$bam_dir"9_20.splitters.merged.bam",\
		$bam_dir"9_5_1.splitters.merged.bam",\
		$bam_dir"9_5_3.splitters.merged.bam",\
		$bam_dir"9_6.splitters.merged.bam",\
	-D $bam_dir"0_1.discordants.merged.bam",\
		$bam_dir"0_2.discordants.merged.bam",\
		$bam_dir"0_2.dis0_4.discordants.merged.bam",\
		$bam_dir"0_2.dis12_3.discordants.merged.bam",\
		$bam_dir"0_2.dis12_4.discordants.merged.bam",\
		$bam_dir"0_2.dis12_5_1.discordants.merged.bam",\
		$bam_dir"0_2.dis1_2.discordants.merged.bam",\
		$bam_dir"0_2.dis13_1.discordants.merged.bam",\
		$bam_dir"0_2.dis13_2.discordants.merged.bam",\
		$bam_dir"0_2.dis13_3.discordants.merged.bam",\
		$bam_dir"0_2.dis13_5_1.discordants.merged.bam",\
		$bam_dir"0_2.dis14_5_1.discordants.merged.bam",\
		$bam_dir"0_2.dis15_5_1.discordants.merged.bam",\
		$bam_dir"0_2.dis16_5_1.discordants.merged.bam",\
		$bam_dir"0_2.dis16_5_2.discordants.merged.bam",\
		$bam_dir"0_2.dis17_5_1.discordants.merged.bam",\
		$bam_dir"0_2.dis17_5_2.discordants.merged.bam",\
		$bam_dir"0_2.dis18_5A.discordants.merged.bam",\
		$bam_dir"0_2.dis19_1B.discordants.merged.bam",\
		$bam_dir"0_2.dis19_3.discordants.merged.bam",\
		$bam_dir"0_2.dis19_4.discordants.merged.bam",\
		$bam_dir"0_2.dis19_5_3.discordants.merged.bam",\
		$bam_dir"0_2.dis19_5.discordants.merged.bam",\
		$bam_dir"0_2.dis20_5.discordants.merged.bam",\
		$bam_dir"0_2.dis21_5_1.discordants.merged.bam",\
		$bam_dir"0_2.dis21_5_2.discordants.merged.bam",\
		$bam_dir"0_2.dis2_1.discordants.merged.bam",\
		$bam_dir"0_2.dis22_5_2B.discordants.merged.bam",\
		$bam_dir"0_2.dis2_5_11.discordants.merged.bam",\
		$bam_dir"0_2.dis2_5_9.discordants.merged.bam",\
		$bam_dir"0_2.dis3_5_15.discordants.merged.bam",\
		$bam_dir"0_2.dis3_5_1.discordants.merged.bam",\
		$bam_dir"0_2.dis3_5_2.discordants.merged.bam",\
		$bam_dir"0_2.dis3_6.discordants.merged.bam",\
		$bam_dir"0_2.dis6_2.discordants.merged.bam",\
		$bam_dir"0_2.dis6_3.discordants.merged.bam",\
		$bam_dir"0_2.dis7_3.discordants.merged.bam",\
		$bam_dir"0_2.dis7_5_4.discordants.merged.bam",\
		$bam_dir"0_2.dis7_5.discordants.merged.bam",\
		$bam_dir"0_2.dis8_5_3.discordants.merged.bam",\
		$bam_dir"0_2.dis9_20.discordants.merged.bam",\
		$bam_dir"0_2.dis9_5_1.discordants.merged.bam",\
		$bam_dir"0_2.dis9_5_3.discordants.merged.bam",\
		$bam_dir"0_2.dis9_6.discordants.merged.bam",\
	-o $lumExp_dir"CNV_multi_sample.vcf"\
	-T $temp_dir







		
	

  	








