#!bin/bash

# Change the header of the fasta file with the sacffold and contigs name
sed -e "s/LRGB01........ Daphnia magna strain Xinb3 //g" GCA_001632505.1_daphmag2.4_genomic.fna | sed 's/, whole genome shotgun sequence//' > Dmagna_scaffID.fna

# Create a file with all the fasta headers
grep '>' GCA_001632505.1_daphmag2.4_genomic.fna > headers.fasta

# Create tab-separaed file with ID anc scaffolds or contigs
sed -e "s/ Daphnia magna strain Xinb3 /\t/g" headers.fasta | sed 's/, whole genome shotgun sequence//'| sed -e "s/^>//" > ID_contigs.txt

# invert the order of the columns
awk '{print $2,$1}' ID_contigs.txt > pattern.txt 

# Substitute all the contigs name to ncbID according to pattern.txt. Remember to make a copy of the file because "-i" change the input file 
while IFS=" " read -r contig ncbID
do
#	sed -i "s/$ncbID/$contig/g" Dmagna_scaffID.gff
	sed -i "s/$ncbID/$contig/g" Dmagna_scaffID.fna.fai 

	echo "from $ncbID to $contig"

done < "pattern.txt"

