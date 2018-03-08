# <i>Daphnia magna</i> Variant Call Analysis
Variant call analysis was performed using the <i>D. magna</i> genome GCA_001632505.1 and the raw sequences of the <i>D. magna</i> samples listed below. Data analysis was performed using three bash scripts for: (i) quality control and alignment, (ii) read group assignment and bam file merging, and (iii) post-mapping analysis and variant calling.  

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
* ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10 (TruSeq3-PE-2.fa is available at: https://github.com/timflutre/trimmomatic/blob/master/adapters/TruSeq3-PE-2.fa)
* LEADING:30 
* TRAILING:30 
* MINLEN:50 

Filtered data was aligned using BWA-mem (Li, 2013) piped to SAMtools view (Li, 2009) and SAMtools sort (Li, 2009) to create a sorted bam file. <i>D. magna</i> genome GCA_001632505.1 was downloaded from NCBI and used as reference.

## Script 2: Read Group Assignments and File Merging
Read group was assigned to each replicate using the Picard tool (http://broadinstitute.github.io/picard/) AddOrReplaceReadGroups as indicated in the raw file names. Then, the three replicates were merged using Picard MergeSamFiles.

## Script 3: Post-mapping Analysis and Variant Calling 
A reference dictionary was created using the Picard tool CreateSequenceDictionary and SortSam was used on merged bam files to sort the reads and create new indexes. Additionally, removal of PCR duplicate generated during library constructions was performed using the Picard tool MarkDuplicate.
Before calling variants using SAMtools mpileup (Li, 2013), local realignment of reads around indels was carried out with the Genome Analysis ToolKit (McKenna et al.,2010) RealignerTargetCreator and IndelRealigner. Ultimately, variants were filtered using BCFtools filter (https://samtools.github.io/bcftools/)


*SortSam settings*
* SORT_ORDER=coordinate
* CREATE_INDEX=true
* VALIDATION_STRINGENCY=SILENT 	

*MarkDuplicate settings*
* ASSUME_SORT_ORDER=coordinate
* REMOVE_SEQUENCING_DUPLICATES=true
* VALIDATION_STRINGENCY=SILENT 

*BCFtools filter settings*
* -i 'MIN(DP)>10 && QUAL>30 && AVG(GQ)>50'


## References
* Bolger, A. M., Lohse, M., & Usadel, B. (2014). Trimmomatic: a flexible trimmer for Illumina sequence data. <i>Bioinformatics</i>, 30(15), 2114-2120.
* Li, H., Handsaker, B., Wysoker, A., Fennell, T., Ruan, J., Homer, N., ... & Durbin, R. (2009). The sequence alignment/map format and SAMtools. <i>Bioinformatics</i>, 25(16), 2078-2079.
* Li, H. (2011). A statistical framework for SNP calling, mutation discovery, association mapping and population genetical parameter estimation from sequencing data. <i>Bioinformatics</i>, 27(21), 2987-2993.
* Li, H. (2013). Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM. <i>arXiv preprint arXiv</i>:1303.3997.
* McKenna, A., Hanna, M., Banks, E., Sivachenko, A., Cibulskis, K., Kernytsky, A., ... & DePristo, M. A. (2010). The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data. <i>Genome research</i>, 20(9), 1297-1303.




