#!/usr/bin/env python

import sys
import argparse
import re

parser = argparse.ArgumentParser(description='Add missing nucleotides to the selected decay reads upto certail length')
parser.add_argument("-i","--input-file", help = "nucleotide sequence file in txt format")
parser.add_argument("-u","--up-range", type = int, default = 100, help = "No of nts upstream of read start. Default: 100")
parser.add_argument("-d","--down-range", type = int, default = 100, help = "No of nts downstream of read start. Default: 100")
parser.add_argument("-o","--output-file", help = "output file in txt format with required nucleotides added")

args = parser.parse_args()

output_file = open(args.output_file,"w+")

with open(args.input_file) as in_file:
    for line in in_file:
        new_line=re.split(" +|\t", line)
        read=(new_line[0])
        start = int(new_line[3])
        upstream = start - args.up_range
        end = int(new_line[4])
        downstream = start + args.down_range
        sequence = (new_line[5].strip())
        sequence = sequence.upper()

        if upstream >= 1 and downstream > end:   # with downstream ends below 3p-ends
            sequence = sequence + str((downstream - end)*"A")
            output_file.write(str(sequence)+ "\n")

        elif upstream < 1 and downstream <= end:   # with upstrean above 5p-ends
            sequence =  str((abs(upstream)+ 1)*"N") + sequence  #1 added to upstream as it is zero based coordinates
            output_file.write(str(sequence)+ "\n")

        elif upstream < 1 and downstream > end:   # odd files where upstream is above start and downstreanm is below end
            sequence = ((abs(upstream))*"N")+ sequence + str((downstream - (end-1))*"A") #end-1 is because of 0-based coordinates in files
            output_file.write(str(sequence)+ "\n")

        else: #For files with 200 ntd within transcript
            if len(sequence)==200:
                output_file.write(str(sequence)+ "\n")
            else:
                print("ERROR in", new_line)

output_file.close()
