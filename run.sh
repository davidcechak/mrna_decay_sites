#!/bin/bash

# Convert genomic coordinates to transcriptomic coordinates
# mkdir -p per_RBP_transcr_BED
# for i in per_RBP_BED/*; do
# 	./dev/genomic-to-transcr-coords/src/bed-genomic-to-transcr-coords.py \
# 		-g ./hg38/genes.gtf \
# 		-b $i \
# 		> per_RBP_transcr_BED/"$(basename -- $i)";
# done



# Convert genomic coordinates to transcriptomic coordinates WITH ALT/REF columns
mkdir -p per_RBP_transcr_BED
for i in per_RBP_BED/*; do
    echo $i
# HERE is the only diff
	./dev/genomic-to-transcr-coords/src/bed-genomic-to-transcr-coords-alt-ref.py \
		-g ./hg38/genes.gtf \
		-b $i \
		> per_RBP_transcr_BED/"$(basename -- $i)";
done
