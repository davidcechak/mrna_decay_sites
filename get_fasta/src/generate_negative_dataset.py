#!/usr/bin/env python

import sys
import argparse
import re
import random

parser = argparse.ArgumentParser()
parser.add_argument("-i","--input-fasta.fai_file", help = "input file of transcript reads in fasta.fai file in bed format")
parser.add_argument("-n","--number_sequences", type = int, default = 1300000, help = "No of nts upstream of read start. Default: 100")
parser.add_argument("-o","--output_random_file", help = "output file in bed format with random coordiantes")
#parser.add_argument("-p","--output_random_file_ii", help = "output file in txt format with required nucleotides added")

args = parser.parse_args()

needed_sequence=args.number_sequences
counts=0

random_file = open(args.output_random_file,"w+")
#random_file_case_2 = open(args.output_random_file_ii,"w+")

print("CASE-I: GENERTAING READS FROM TRANSCRIPTS RANDOMLY FROM START TO END")

while (counts <= needed_sequence):
    with open(args.input_fasta_file) as transcripts_file:
        for transcript in transcripts_file:
            len_seq=len(transcript)
            random_start=random.randint(1,len_seq+1)
            new_line=re.split("\t", transcript)
            read=(new_line[0])
            start = int(new_line[3])
            len_seq=len(transcript)
            # had a debugging error IndentationError: unexpected indent
                # random_start=random.randint(1,len_seq+1)
                # rand_sequence=transcript[random_start:start_plus_200].upper()
                # rand_seq_final=rand_sequence.strip()
                # random_file.write(str(rand_seq_final)+ "\n")
                # counts+=1
            random_start=random.randint(1,len_seq+1)
            rand_sequence=transcript[random_start:start_plus_200].upper()
            rand_seq_final=rand_sequence.strip()
            random_file.write(str(rand_seq_final)+ "\n")
            counts+=1


# print("CASE-II: GENERTAING READS FROM TRANSCRIPTS WITH NUCLEOTIDES GREATER THAN 300 & LAST 200 NUCLEOTIDES IGNORED AS START")
#
# while (counts <= needed_sequence):
#     with open(args.input_fasta_file) as transcripts_file:
#         for transcript in transcripts_file:
#             if ">" in transcript:
#                 pass
#             else:
#                 len_seq=len(transcript)
#                 if len_seq > 300:
#                     random_start=random.randint(1,len_seq-201)
#                     start_plus_200=random_start+200
#                     rand_sequence=transcript[random_start:start_plus_200].upper()
#                     rand_seq_final=rand_sequence.strip()
#                     output_file_case_2.write(str(rand_seq_final)+ "\n")
#                     counts+=1
#                 else:
#                     pass
# #print("Total : ", counts)
