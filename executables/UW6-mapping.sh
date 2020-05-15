#! /bin/bash 

BOWTIE=/opt/bifxapps/bowtie2-2.2.1
SAMTOOLS=/opt/bifxapps/samtools-1.9/bin

ref=$1
metagenome=$2
out=$3

$BOWTIE/bowtie2 -p 4 -x $ref -q $metagenome > $out

# Get soft reads
$SAMTOOLS/samtools view -f 2 -F 256 -bS $out > ${out%.sam}.soft.bam

# Sort soft reads
$SAMTOOLS/samtools sort ${out%.sam}.soft.bam -o ${out%.sam}.soft.sorted.bam

# Get hard reads
$SAMTOOLS/samtools view -bS -f 2 -F 3848 $out > ${out%.sam}.hard.bam

# Sort hard reads 
$SAMTOOLS/samtools sort ${out%.sam}.hard.bam -o ${out%.sam}.hard.sorted.bam