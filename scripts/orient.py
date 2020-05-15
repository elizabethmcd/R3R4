#!/usr/bin/python3.6
# Alec Bahcheli
# name is dnaA or whatever the dnaa locus is called in the gff file that you found it as

import sys, re, math
dnaa = sys.argv[1]
finalassembly = sys.argv[2]
output = sys.argv[3]
name = sys.argv[4]
shiftvalue = open(sys.argv[5], "w")

dnaa=open(dnaa).readlines()

genome = open(finalassembly)
# open a file for outputting the oriented genome
orientedgenome = open(output, "w")
# spots will be the potential ori-starts
spot = []
# define the raw genome
raw = ""
# get the dnaa gene loci
for line in dnaa:
    if name in line:
        line = line.split("\t")
        spot.append(int(line[4]) - 1)

if len(spot) != 0:
    for line in genome:
        if ">" in line:
            line = re.sub("circular=true_pilon",'',line)
            orientedgenome.write(line.strip("\n") + "\n")
        else:
            raw = raw + line.strip("\n")
    # there are usually 60nts per line in fasta sequences
    numlines = math.ceil(len(raw) / 60)
    result = raw[int(spot[0]):] + raw[:int(spot[0])]
    for i in range(numlines):
        orientedgenome.write(result[(i * 60): min(((i + 1) * 60), len(raw))] + "\n")
    orientedgenome.close()
else:
    for line in genome:
        if ">" in line:
            line = re.sub("circular=true_pilon",'',line)
            orientedgenome.write(line.strip("\n") + "\n")
        else:
            raw = raw + line.strip("\n")
    # there are usually 60nts per line in fasta sequences
    numlines = math.ceil(len(raw) / 60)
    result = raw[int(spot[0]):] + raw[:int(spot[0])]
    for i in range(numlines):
        orientedgenome.write(result[(i * 60): min(((i + 1) * 60), len(raw))] + "\n")
    orientedgenome.close()
genome.close()
shiftvalue.write(str(spot[0]) + "\t" + str(len(raw)) + "\n")
shiftvalue.close()