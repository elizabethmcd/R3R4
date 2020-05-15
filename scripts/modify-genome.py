#!/usr/bin/python3.6
# Alec Bahcheli
# Rebuild the genome to include ~80kb from the end copied to the front

import sys, math

finaloriented = sys.argv[1]
output = sys.argv[2]
# get oriented genome
genome = open(finaloriented).readlines()
# prepare reversed genome
output = open(output, "w")
seq = ""
for line in genome:
    if ">" in line:
        line = line.strip("\n")
        line = line + "-reversed" + "\n"
        output.write(line)
    else:
        seq = seq + line.strip("\n")
# determine where to cut genome (half of genome to lowest thousand)
cut = math.floor((math.ceil(len(seq) / 1000)) / 2)
start = cut * 1000
# re-write the genome from the half-way point for mapping
result = seq[start:] + seq[:start]
# fasta files usually have 60nts per line
numlines = math.ceil(len(result) / 60)
for i in range(numlines):
    output.write(result[(i * 60): min(((i + 1) * 60), len(result))] + "\n")
output.close()