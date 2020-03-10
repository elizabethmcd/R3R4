#! /usr/bin/bash

# Queue metatranscriptomic mapping to reference genomes/MAGs with a metadata file specifying the sample name and filenames of the R1 and R2 files
# Metadata file formatted as FilenameR1 FilenameR2 Folder/Sample Name Such as:
# SRX4072504.qced.R1.fastq	SRX4072504.qced.R2.fastq	20120800_P2M
# in a tab delimited file
# usage: bash queue-kallisto.sh $index $MetadataFile
# Make sure you have indexed your reference genomes first before performming mapping, as you need to provide the resulting index file from "kallisto index"

index="$1"
MetadataFile="$2"

echo "Running kallisto for $MetadataFile !"

while read -r a b c; do
	echo $c $a $b;
	kallisto quant -i "$index" -o "$c" "$a" "$b";
done < "$MetadataFile"

echo "Finished running kallisto for $MetadataFile !"
