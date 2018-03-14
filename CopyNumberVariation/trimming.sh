#!/bin/sh

module load apps/trimmomatic/v0.32
module load apps/fastqc/v0.11.2

### Set up working directories for trimmomatic
HomeDir='/rds/projects/2017/orsinil-bioinfoproj/CNV/'
mkdir -p $HomeDir'trimReads'
mkdir -p $HomeDir'unpairReads'
trimReads=$HomeDir'trimReads/'
unpairReads=$HomeDir'unpairReads/'
samples=$HomeDir'RawFwdFqLists.txt'

### Run trimmomatic
while IFS= read -r sample ; do

	file1=$sample
	fileName1=$(basename "$file1")
	
	file2=${sample//R1.fastq.gz/R2.fastq.gz}
	fileName2=$(basename "$file2")


	trim1=${fileName1//"-R1.fastq.gz"/"_1_trim.fastq.gz"}
	trim2=${fileName2//"-R2.fastq.gz"/"_2_trim.fastq.gz"}

	unpair1=${fileName1//".R1.fastq.gz"/"_1_trim.fastq.gz"}
	unpair2=${fileName2//".R2.fastq.gz"/"_2_trim.fastq.gz"}

		echo $sample
		echo $file1
		echo $fileName1
		echo $file2
		echo $fileName2
		echo $trim1
		echo $trim2
		echo $unpair1
		echo $unpair2
		echo ""
		 
	#java -Xms10g -Xmx20G -jar  ${TMATIC_ROOT}/trimmomatic-0.32.jar PE $file1 $file2 $trimReads$trim1 $unpairReads$unpair1 $trimReads$trim2 $unpairReads$unpair2 ILLUMINACLIP:adapter_seqs.fa:2:30:10 LEADING:30 TRAILING:30 MINLEN:50 -threads 4

	trimmomatic PE $file1 $file2 $trimReads$trim1 $unpairReads$unpair1 $trimReads$trim2 $unpairReads$unpair2 ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10 LEADING:30 TRAILING:30 MINLEN:50


done < "$samples"



### Set up working directories for fastqc
mkdir -p $HomeDir'trimFQC'
trimFQC=$HomeDir'trimFQC/'
ls $trimReads*.fastq.gz > $HomeDir'TrimFqLists.txt'
trimFastqs=$HomeDir'TrimFqLists.txt'

### Run fastqc
while IFS= read -r trimFastqs ; do
	
	echo $trimFasrqs
	
	fastqc $trimFastqs --threads 4 --outdir $trimFQC

done < "$trimFastqs"




