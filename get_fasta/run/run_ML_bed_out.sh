#!/bin/bash

FILTER_CODING_REGIONS=0
FILTER_PROMINENT_REGIONS=0
FILTER_ADAPTIVE_PROMINENT=1
RECOMPUTE_TRAIN_TEST_BAMS=0
region_diameter=5
threshold=3

# version of adaptive filtering
# version="cds_local.max.ds2"
version="adaptive"
# version="cds_local.max.ds2to11"


# location (folder) of the result
result_folder='/original'
if test $FILTER_PROMINENT_REGIONS = 1 && test $FILTER_CODING_REGIONS = 1; then
    result_folder='/coding_n_not_prominent'
elif test $FILTER_PROMINENT_REGIONS = 1; then
    result_folder='/not_prominent'        
elif test $FILTER_CODING_REGIONS = 1; then
    result_folder='/coding'
elif test $FILTER_ADAPTIVE_PROMINENT = 1; then
    result_folder="/${version}/adaptive"
else
    result_folder='/original'
fi
echo $result_folder
echo 'region_diameter = '$region_diameter 
echo 'threshold = '$threshold



samples_path="../samples/"
data_path="../data"
# Input file containing the reads
reads_file_name="reads.1.sanitize.noribo.toTranscriptome.sorted.bam"
# reads_file_name="reads.1.sanitize.noribo.toTranscriptome.sorted.uxi.bam"

#DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
#SAMPLE_DIR="$results"

samples=(
	# "hsa.dRNASeq.HeLa.polyA.REL5.long.1"
    # "hsa.dRNASeq.HeLa.polyA.REL5.long.2"
	# "hsa.dRNASeq.HeLa.polyA.REL5.long.3"
	# "hsa.dRNASeq.HeLa.polyA.REL5.long.4"
	# "hsa.dRNASeq.HeLa.polyA.REL5.long.5"
	# "hsa.dRNASeq.HeLa.polyA.REL5.long.6"
	# "hsa.dRNASeq.HeLa.polyA.REL5.long.7"
	# "hsa.dRNASeq.HeLa.polyA.REL5.long.8"
	# "hsa.dRNASeq.HeLa.polyA.REL5.long.9"
	# "hsa.dRNASeq.HeLa.polyA.REL5.long.10"
	# "hsa.dRNASeq.HeLa.polyA.REL5.long.11"
    # "ds2to11"
    # "ds1ds2"
    "hsa_dRNA_HeLa_NoArs_polyA_5P_1"

)

for name in "${samples[@]}"; do
    mkdir -p results/$name/$result_folder/
done


for name in "${samples[@]}"; do
    if test $RECOMPUTE_TRAIN_TEST_BAMS = 1; then
        ### Filter chromosome 1 sequences as a test set
        # get transcript identifications belonging to chromosome 1 from the GTF (Gene transfer format) annotation file
        awk '$1 == "chr1" && $3 == "exon" {for (i=1; i<=NF; i++) if ($i=="transcript_id") {print $(i+1)} }' ../data/annotation.gtf | tr -d '";' | sort | uniq > chr1_transcripts.txt

        chmod 777 chr1_transcripts.txt
        # Convert the BAM file to a SAM file
        samtools view -h -o reads.sam $samples_path$name/$reads_file_name

        # filter the reads that belong to chr1
        # adjusts the identifiers in chr1_transcripts.txt to match the format in the reads.sam file
        grep -F -f <(sed 's/\..*//' chr1_transcripts.txt) reads.sam > reads_chr1.sam

        # convert the SAM file back to a BAM file for further processing
        samtools view -bS reads_chr1.sam > $samples_path$name"/"reads_chr1.bam

        reads_file_name_test="reads_chr1.bam"
        # create an index file for it 
        samtools index $samples_path$name"/"$reads_file_name_test

        # TODO remove intermediate file when I do not need them anymore
        # rm reads.sam reads_chr1.sam

        # filter the reads that DO NOT belong to chr1
        grep -v -F -f <(sed 's/\..*//' chr1_transcripts.txt) reads.sam > reads_no_chr1.sam

        # convert the SAM file back to a BAM file for further processing
        samtools view -bS reads_no_chr1.sam > $samples_path$name"/"reads_no_chr1.bam

        reads_file_name_train="reads_no_chr1.bam"
        # create an index file for it 
        samtools index $samples_path$name"/"$reads_file_name_train
        # TODO remove intermediate file when I do not need them anymore
        # rm reads.sam reads_no_chr1.sam
    else
        reads_file_name_test="reads_chr1.bam"
        reads_file_name_train="reads_no_chr1.bam"
    fi
done


for name in "${samples[@]}"; do
    if test $FILTER_CODING_REGIONS = 1; then
        echo ">>> FILTER: CODING REGIONS <<<" 
        # >>> FILTER: LEAVE ONLY READS WHOSE STARTING NUCLEOTIDE IS LOCATED WITHIN THE CDS REGIONS <<<"
        for split in "train" "test"; do
            reads_file_name_var="reads_file_name_$split"
            previous_reads_file_name=${!reads_file_name_var}
            result_name="coding_filtered_${split}.bam"
            
            samtools index $samples_path$name"/"$previous_reads_file_name
            # Filter the BAM file based on CDS regions
            bedtools intersect \
                -abam $samples_path$name"/"$previous_reads_file_name \
                -b run/genic_elements.cds.bed \
                -wa \
                -u \
                > $samples_path$name"/"$result_name
            # Index the BAM which contains only CDS 
            samtools index $samples_path$name"/"$result_name
            echo $samples_path$name"/"$previous_reads_file_name

            eval reads_file_name_$split=$result_name
            
            # Count the number of reads in the original BAM file
            echo "original # of reads:"
            samtools view -c $samples_path$name"/"$previous_reads_file_name

            # Count the number of reads in the filtered BAM file
            echo "filtered # of reads:"
            samtools view -c $samples_path$name"/"${!reads_file_name_var}
            
            echo $reads_file_name_var
            echo $samples_path$name/${!reads_file_name_var}
        done
    fi
done

# 
threshold_ratio=0.01  # define an appropriate ratio for setting the threshold

for name in "${samples[@]}"; do
    if test $FILTER_PROMINENT_REGIONS = 1; then
        echo ">>> FILTER: PROMINENT REGIONS <<<" 
        # >>> FILTER: LEAVE ONLY READS WHOSE STARTING NUCLEOTIDE IS LOCATED IN A DENSE AREA <<<
        for split in "train" "test"; do
            reads_file_name_var="reads_file_name_$split"
            previous_reads_file_name=${!reads_file_name_var}
            result_name="prominent_filtered_${split}.bam"
            # Convert the BAM file to a BED file
            bedtools bamtobed \
                -i $samples_path$name/$previous_reads_file_name > $samples_path$name/reads_$split.bed
            wait

            # Extract the starting positions of the reads and extend the regions by X/2 on both sides
            awk -v X=$region_diameter '{start = ($2 - int(X/2) >= 0) ? $2 - int(X/2) : 0; end = $2 + int(X/2) + 1; print $1 "\t" start "\t" end}' $samples_path$name/reads_$split.bed > $samples_path$name/read_starts_extended_$split.bed
            wait

            # Calculate the read start density for the extended regions
            bedtools intersect \
                -a $samples_path$name/read_starts_extended_$split.bed \
                -b $samples_path$name/reads_$split.bed \
                -c > $samples_path$name/read_start_density_extended_$split.bed
            wait

            # Filter the regions based on the density threshold
            awk -v threshold=$threshold '{if ($4 >= threshold) print $0}' $samples_path$name/read_start_density_extended_$split.bed > $samples_path$name/dense_regions_$split.bed

            # Filter the reads based on the dense regions
            bedtools intersect -a $samples_path$name/$previous_reads_file_name -b $samples_path$name/dense_regions_$split.bed \
                -u > $samples_path$name/$result_name
            
            eval reads_file_name_$split=$result_name
            
            # Create index file
            samtools index $samples_path$name"/"${!reads_file_name_var}

            # Count the number of reads in the original BAM file
            echo "original # of reads:"
            samtools view -c $samples_path$name"/"$previous_reads_file_name

            # Count the number of reads in the filtered BAM file
            echo "filtered # of reads:"
            samtools view -c $samples_path$name/${!reads_file_name_var}
            wait
            
            echo $samples_path$name/${!reads_file_name_var}
        done
    fi
done


if test $FILTER_ADAPTIVE_PROMINENT = 1; then
    echo "FILTER_ADAPTIVE_PROMINENT"
    reads_file_name_test="adaptive_filtered_reads_test.sorted.bam"
    reads_file_name_test="/${version}/adaptive_filtered_reads_test.sorted.bam"
    
    reads_file_name_train="adaptive_filtered_reads_train.sorted.bam"
    reads_file_name_train="/${version}/adaptive_filtered_reads_train.sorted.bam"
fi


echo ">>> EXTRACTING THE UPSTREAM AND DOWNSTREAM NUCLEOTIDE COORDINATES BASED ON mRNA DECAY PROFILE <<<"
for name in "${samples[@]}"; do
    for split in "train" "test"; do
        reads_file_name_var="reads_file_name_$split"
            
        mkdir -p results/$name$result_folder/trans_coords
        echo ">>> EXTRACTING TRANSCRIPT NUCLEOTIDE COORDINATES IN MODIFIED BED FILES <<<"
        echo ">>>  Working for" $name
        src/extract_5p_decay_mrna_coordinates_for_ML.py \
            -b $samples_path$name${!reads_file_name_var} \
            -x results/$name$result_folder/trans_coords/upstream-5P-end-transcripts_$split.bed \
            -y results/$name$result_folder/trans_coords/within_transcripts_$split.bed \
            -z results/$name$result_folder/trans_coords/downstream-3P-end-transcripts_$split.bed \
            -w results/$name$result_folder/trans_coords/odd_transcripts_$split.bed \
            -v results/$name$result_folder/trans_coords/extra_transcript_$split.bed

        echo ">>> COLLECTING COORDINATES VIA BEDTOOLS<<<"
        bedtools bamtobed \
        -i $samples_path$name${!reads_file_name_var} \
         > results/$name$result_folder/trans_coords/bamtobed_all_file_$split.bed
     done
done


echo ">>> COLLECTING NUCLEOTIDE SEQUENCE FROM BED-COORDINATES <<<"
for name in "${samples[@]}"; do
    for split in "train" "test"; do
        mkdir -p results/$name$result_folder/decay_seqs_fasta
        mkdir -p results/$name$result_folder/decay_seqs_txt
        mkdir -p results/$name$result_folder/decay_seqs__200ntd_txt
        mkdir -p results/$name$result_folder/random_neg_dataset

        echo ">>>  WORKING FOR" $name "<<<"
        echo ">>> COLLECTING FASTA SEQUENCES FOR COORDINATES THAT LIE WITHIN TRANSCRIPTS <<<"
        #bedtools getfasta -fi test.fa -bed new.bed -fo outfile.fa
        bedtools getfasta \
            -fi $data_path/hg38/transcripts.fa \
            -bed results/$name$result_folder/trans_coords/within_transcripts_$split.bed \
            -fo results/$name$result_folder/decay_seqs_fasta/within_transcripts_fasta_$split.fa
            # THIS is the data used for training later -- decay_seqs_fasta/within_transcripts_fasta.fa

            grep "^[[:alpha:]]" results/$name$result_folder/decay_seqs_fasta/within_transcripts_fasta_$split.fa \
                > results/$name$result_folder/decay_seqs_txt/within_transcripts_ntd_seqs_$split.txt

        echo ">>> COLLECTING FASTA SEQUENCES FOR COORDINATES THAT LIE UPSTREAM 5P-END OF TRANSCRIPTS <<<"
        bedtools getfasta \
            -fi $data_path/hg38/transcripts.fa \
            -bed results/$name$result_folder/trans_coords/upstream-5P-end-transcripts_$split.bed \
            -fo results/$name$result_folder/decay_seqs_fasta/upstream-5P-end-transcripts_fasta_$split.fa
            # THIS is the data used for training later -- decay_seqs_fasta/upstream-5P-end-transcripts_fasta.fa

        grep "^[[:alpha:]]" results/$name$result_folder/decay_seqs_fasta/upstream-5P-end-transcripts_fasta_$split.fa \
            > results/$name$result_folder/decay_seqs_txt/upstream-5P-end-transcript_ntd_seqs_$split.txt

        echo ">>> COLLECTING FASTA SEQUENCES FOR COORDINATES THAT LIE DOWNSTREAM 3P-END OF TRANSCRIPTS <<<"
        bedtools getfasta \
            -fi $data_path/hg38/transcripts.fa \
            -bed results/$name$result_folder/trans_coords/downstream-3P-end-transcripts_$split.bed \
            -fo results/$name$result_folder/decay_seqs_fasta/downstream-3P-end-transcripts_fasta_$split.fa
            # THIS is the data used for training later -- decay_seqs_fasta/downstream-3P-end-transcripts_fasta.fa

            grep "^[[:alpha:]]" results/$name$result_folder/decay_seqs_fasta/downstream-3P-end-transcripts_fasta_$split.fa \
                > results/$name$result_folder/decay_seqs_txt/downstream-3P-end-transcripts_fasta_ntd_seqs_$split.txt
    done
done

echo ">>> FIANLISED ALL RUNS <<<"
