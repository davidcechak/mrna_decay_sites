#!/usr/bin/env python

import argparse
import re
from intervaltree import Interval, IntervalTree

class Region:

    def __init__(self, chrom, strand, start, end):
        self.chrom = chrom
        self.strand = strand
        self.start = start
        self.end = end

    def len(self):
        return self.end - self.start

class Transcript:

    def __init__(self, transcript_id):
        self.id = transcript_id
        self.exons = []
        self.start = None
        self.end = None
        self.chrom = None
        self.strand = None
        self.cds = Region(None, None, None, None)
        self.gid = None

    def add_exon(self, region):
        self.exons.append(region)
        if self.start == None or region.start < self.start:
            self.start = region.start
        if self.end == None or region.end > self.end:
            self.end = region.end
        if self.chrom == None:
            self.chrom = region.chrom
        if self.strand == None:
            self.strand = region.strand

    def add_terminal_codon(self, region):
        if self.cds.start == None or region.start < self.cds.start:
            self.cds.start = region.start
        if self.cds.end == None or region.end > self.cds.end:
            self.cds.end = region.end

    def is_region_in_exons(self, chrom, strand, start, end):
        if self.chrom != chrom:
            return False
        if self.strand != strand:
            return False
        for e in self.exons:
            if start >= e.start and end <= e.end:
                return True
        return False

def gtf_to_transcript_list(gtf):
    transcript_id_map = {}
    with open(gtf) as infile :
        for line in infile :
            if line[0] == "#":
                continue
            cols = line.strip().split('\t')
            chrom = cols[0]
            line_type = cols[2]
            start = int(cols[3]) - 1 # GTF is 1-based; convert to 0-based [start, end)
            end = int(cols[4])
            strand = cols[6] #strand is + or -
            rest = cols[8]
            gid = re.findall(r"gene_id\s+\"(\w+)\"", rest)[0]
            tid = re.findall(r"transcript_id\s+\"(\w+)\"", rest)[0]

            # Find the correspoding transcript for this line.
            tr = None
            if tid in transcript_id_map:
                tr = transcript_id_map[tid]
            else:
                tr = Transcript(tid)
                tr.gid = gid
                transcript_id_map[tid] = tr

            region = Region(chrom, strand, start, end)
            if line_type == "exon":
                tr.add_exon(region)
            elif line_type == "start_codon" or line_type == "stop_codon":
                tr.add_terminal_codon(region)

    transcripts = transcript_id_map.values()
    return transcripts


parser = argparse.ArgumentParser()
parser.add_argument("-b","--BED-file", help = "BED with genomic coords")
parser.add_argument("-g","--GTF-file", help = "GTF with gene models")
args = parser.parse_args()

# Read the genomic coordinates of the transcripts from the GTF and create
# transcript objects.
transcripts = gtf_to_transcript_list(args.GTF_file)

# Add the transcripts in a interval tree for efficient range queries.
tree = IntervalTree()
for t in transcripts:
    interv = Interval(t.start, t.end, t)
    tree.add(interv)

#Open and start reading the BED file.
with open(args.BED_file) as infile:
    for line in infile :
        #(chrom, start, end, name, score, strand) = line.strip().split('\t')[0:6]
        # (chrom, start, end, name, score, strand, ref_alt) = line.strip().split('\t')[0:7]
        (chrom, start, end, name, score, strand, ref_alt, info) = line.strip().split('\t')[0:8]
        start = int(start)
        end = int(end)
        chrom = chrom.replace('chr','')

        overlapping_intervals = tree[start:end]
        for interval in overlapping_intervals:
            t = interval.data # the transcript is stored in the interval data

            if not t.is_region_in_exons(chrom, strand, start, end):
                continue

            converted_start = 0
            converted_end = 0
            for e in t.exons:
                if strand == "+":
                    if start >= e.end :
                        converted_start += e.len()
                    elif start > e.start and start < e.end:
                        converted_start += start - e.start
                    if end >= e.end:
                         converted_end += e.len()
                    elif end < e.end and end > e.start:
                         converted_end += end - e.start
                elif strand == "-":
                    if start <= e.start:
                         converted_end += e.len()
                    elif start > e.start and start < e.end:
                         converted_end += e.end - start
                    if end <= e.start:
                         converted_start += e.len()
                    elif end > e.start and end < e.end:
                         converted_start += e.end - end

            #print("\t".join([t.id, str(converted_start), str(converted_end), name, score, "+"]))
            print("\t".join([t.id, str(converted_start), str(converted_end), name, score, "+", ref_alt, info]))

