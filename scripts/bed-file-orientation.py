#!/usr/bin/python3.6
# Alec Bahcheli

import sys, re, math
shift_size = open(sys.argv[1]).read()
trna = sys.argv[2]
rrna = sys.argv[3]
cdspos = sys.argv[4]
cdsneg = sys.argv[5]

separate = shift_size.split("\t")
shift = int(separate[0])
size = int(separate[1])
list1=[trna, rrna, cdspos, cdsneg]

for x in list1:
    file1 = open(x).readlines()
    file1_list = []
    y = 0
    for line in file1:
        if y > 0:
            original = line.split("\t")
            original[1] = str(int(original[1]) - shift)
            if int(original[1]) < 0:
                original[1] = str(int(original[1]) + size)
            original[2] = str(int(original[2]) - shift)
            if int(original[2]) < 0:
                original[2] = str(int(original[2]) + size)
            file1_list.append("\t".join(original))
            y += 1
        else:
            file1_list.append(line) 
            y += 1
    file1_results = open(x, "w")
    for spot in file1_list:
        file1_results.write(spot)
    file1_results.close()