#! /bin/bash

###################
# Filtering metagenomic samples with BBduk part of the BBtools package
# For use on WEI GLBRC servers running HT Condor
# Elizabeth McDaniel 
##################

# set path where fastp is installed in local home directory bin
BBPATH=/opt/bifxapps/bbmap-38.32/

# queueing r1 r2 metagenomic reads and output folder/file names
in=$1
out=$2

$BBPATH/bbduk.sh in=$in out=$out qtrim=r trimq=10 maq=10
