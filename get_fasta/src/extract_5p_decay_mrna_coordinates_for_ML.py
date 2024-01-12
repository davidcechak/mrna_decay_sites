#!/usr/bin/env python

# python programme to exract the up-stream and downstream coordinated
#of mRNA decay fragments w.r.t. 5-prime reads

import pysam
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-u","--up-range", type = int, default = 100,
            help = "No of nts upstream of read start. Default: 100")
parser.add_argument("-d","--down-range", type = int, default = 100,
            help = "No of nts downstream of read start. Default: 100")
parser.add_argument("-b","--bam-file",
            help = "RNAseq file in bam format")
parser.add_argument("-x","--upstream-5P-end-transcripts",
            help = "output file in bed format for upstream-5P-end-transcripts")
parser.add_argument("-y","--within_transcripts",
            help = "output file in bed format within_transcripts")
parser.add_argument("-z","--downstream-3P-end-transcripts",
            help = "output file in bed format for downstream-3P-end-transcripts")
parser.add_argument("-w","--odd_transcript",
            help = "output file in bed format for odd-transcripts with cal-start upstream 5'-end and cal-end downstream 3-end")
parser.add_argument("-v","--extra_transcripts",
            help = "output file in bed format extra_transcripts")

args = parser.parse_args()


upstream_5P_ends= open(args.upstream_5P_end_transcripts,"w+")
within_transcripts= open(args.within_transcripts,"w+")
downstream_3P_ends= open(args.downstream_3P_end_transcripts,"w+")
odd_transcript= open(args.odd_transcript,"w+")
extra_transcripts= open(args.extra_transcripts,"w+")


bamfile = pysam.AlignmentFile(args.bam_file, "rb")
#bamf = bamfile.fetch(until_eof=True)
for seq in bamfile:
    if seq.is_unmapped:
        continue

    chrom = seq.reference_name     #chromosome/Transcript
    start = seq.reference_start      #strart
    #end = seq.reference_end - 1 # need to subtract 1 {0 based calculation in pysam}
    end = seq.reference_end

    strand = "+"
    pos = start

    up_stream_end_taken = pos - args.up_range
    down_stream_end_taken = pos + args.down_range

    if up_stream_end_taken < 1 and down_stream_end_taken <= end:
        up_stream_end_taken = 1
        upstream_5P_ends.write(str(chrom) + "\t" + str(up_stream_end_taken) + "\t" + str(down_stream_end_taken) + "\t" + str(pos) + "\t" + str(end) + "\n")

    elif up_stream_end_taken < 1 and down_stream_end_taken >= end:
        up_stream_end_taken = 1
        odd_transcript.write(str(chrom) + "\t" + str(up_stream_end_taken) + "\t" + str(end) + "\t" + str(pos) + "\t" + str(end) + "\n")

    elif up_stream_end_taken > 0 and down_stream_end_taken >= end:
        #down_stream_end_taken = end
        downstream_3P_ends.write(str(chrom) + "\t" + str(up_stream_end_taken) + "\t" + str(end) + "\t" + str(pos) + "\t" + str(end) + "\n")

    elif up_stream_end_taken > 0 and down_stream_end_taken <= end:
        #down_stream_end_taken = end
        within_transcripts.write(str(chrom) + "\t" + str(up_stream_end_taken) + "\t" + str(down_stream_end_taken) + "\t" + str(pos)+  "\t" + str(end)+ "\n")

    else: #checking for other discrepencies in data or processing
         extra_transcripts.write(str(chrom) + "\t" + str(up_stream_end_taken) + "\t" + str(down_stream_end_taken) + "\t" + str(pos)+  "\t" + str(end)+ "\n")

upstream_5P_ends.close()
within_transcripts.close()
downstream_3P_ends.close()
odd_transcript.close()
extra_transcripts.close()

    #print (str(chrom) + "\t" + str(up_stream_end_taken) + "\t" + str(down_stream_end_taken) + "\t" + str(pos) + "\t" + str((down_stream_end_taken)-(up_stream_end_taken)))
