#! /bin/bash

######################
# SPAdes coassembly of metagenomic reads
# Elizabeth McDaniel
######################

# Make sure reads have been filtered previously before assembly - garbage in = garbage out
outdir=/home/GLBRCORG/emcdaniel/EBPR/R3R4/assemblies/R3R4_spades_coassembly

# Python path
#PYTHONPATH=/opt/bifxapps/bin/python3.4
    # use specific 3.4 python path for running SPADES
SPADESPATH=/opt/bifxapps/SPAdes-3.9.0-Linux/bin/

# Run spades
/opt/bifxapps/bin/python3.4 $SPADESPATH/spades.py -t 16 -m 500 -k 21,33,55,77,99,127 --dataset /home/GLBRCORG/emcdaniel/EBPR/R3R4/metadata/R3R4-coassembly.yaml -o $outdir
