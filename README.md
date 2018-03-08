# <i>Daphnia magna</i> Variant Call Analysis
Variant call analysis was performed using the <i>D. magna</i> reference genome GCA_001632505.1 and the raw sequences of the D. magna samples listed below. Data analysis was performed using three bash scripts for: (i) quality control and alignment, (ii) read group assignment and bam file merging, and (iii) post-mapping analysis and variant calling.  

### File list:
* Dmagna1-0_1-L002-R1.fastq.gz
* Dmagna1-0_1-L002-R2.fastq.gz
* Dmagna12-6_2-L002-R1.fastq.gz
* Dmagna12-6_2-L002-R2.fastq.gz
* Dmagna19-9_20-L002-R1.fastq.gz
* Dmagna19-9_20-L002-R2.fastq.gz
* Dmagna31-0_1-L001-R1.fastq.gz
* Dmagna31-0_1-L001-R2.fastq.gz
* Dmagna42-6_2-L001-R1.fastq.gz
* Dmagna42-6_2-L001-R2.fastq.gz
* Dmagna49-9_20-L001-R1.fastq.gz
* Dmagna49-9_20-L001-R2.fastq.gz
* Dmagna61-0_1-L002-R1.fastq.gz
* Dmagna61-0_1-L002-R2.fastq.gz
* Dmagna72-6_2-L002-R1.fastq.gz
* Dmagna72-6_2-L002-R2.fastq.gz
* Dmagna79-9_20-L002-R1.fastq.gz
* Dmagna79-9_20-L002-R2.fastq.gz

## Script 1: Quality Control and Alignment
Raw data quality control was performed using Trimmomatic (Bolger et al., 2014) and FastQC (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

*Trimmomatic settings*
* ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10 
* LEADING:30 
* TRAILING:30 
* MINLEN:50 




using FastQC and MultiQC
Failed tests: Sequence Quality Histograms (11/18); Per Sequence GC Content (18/18); Adapter Content (17/18)

Raw data quality filtering using Trimmomatic
CROP:225 (Keep the first 225 nt)
ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10 (Remove adapters present included in the fasta file TruSeq3-PE-2.fa)
LEADING:10  (Trim 3’ bases with less than 10 quality score)
TRAILING:10 (Trim 5’ bases with less than 10 quality score)
SLIDINGWINDOW:50:20 (Remove reads with an average of quality score lower than 20 over 50 nt) 
MINLEN:36 (Remove reads shorter than 36 nt)

Filtered data quality control using FastQC and MultiQC
Failed tests: Per Sequence GC Content (18/18)

Filtered data alignment against the reference genome using BWA mem piped to samtools view and samtools sort to create a sorted bam file.
The option -M to Mark shorter split hits as secondary should have been added (for Picard compatibility).

## Script 2: Read Group and File Merging
Read group was assigned to each replicate using Picard AddOrReplaceReadGroups
Read group identifier: unique replicate name (e.g. Dmagna1) . 
DNA preparation library identifier: the same as the platform number assuming that all the libraries in the same lane were prepared together 
Platform: illumina 
Platform Unit: flow cell lane number 
Sample: sample name (e.g. Dmagna0_1)

The three replicates were merged using Picard MergeSamFiles 

## Script 3: Post-mapping Analysis and Variant Calling 
A reference dictionary was created using Picard CreateSequenceDictionary 

Picard SortSam was used on merged bam files to sort the reads and create new indexes.
SORT_ORDER=coordinate
CREATE_INDEX=true
VALIDATION_STRINGENCY=SILENT 	
Removal of eventual PCR duplicate reads using Picard MarkDuplicate (Failed for limited system resources of the virtual machine).
ASSUME_SORT_ORDER=coordinate
REMOVE_SEQUENCING_DUPLICATES=true
VALIDATION_STRINGENCY=SILENT 

Define the intervals for local realignment using GATK RealignerTargetCreator on the bam file treated by Picard MarkDuplicate. 
Local realignment of reads around indels defined by GATK RealignerTargetCreator using GATK IndelRealigner on the bam file treated by Picard MarkDuplicate.

Variant calling using Samtools mpileup.


