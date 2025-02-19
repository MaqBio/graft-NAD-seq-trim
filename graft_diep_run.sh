#!/bin/bash

# Setup main dir
BASE_DIR="/Data05/quanma/graft_NADseq"

# go through the whole files with Sample_
for dir in "$BASE_DIR"/Sample_*; do
    # mksure the file exist
    if [[ -d "$dir" ]]; then
        # search *_combined_R2.fastq.gz file
        FASTQ_FILE=$(find "$dir" -type f -name "*_combined_R2.fastq.gz")

        # get fileID（fileID are T1, T2, W1, W2）
        FILE_ID=$(basename "$dir" | awk -F'-' '{print $2}')

        # mksure FASTQ exist
        if [[ -n "$FASTQ_FILE" ]]; then
            echo "Processing $FASTQ_FILE with fileID $FILE_ID"
            cd "$dir" || exit 1
            ./graft-nad-seq.sh "$FASTQ_FILE" "$dir" "$FILE_ID"
            cd "$BASE_DIR" || exit 1
        else
            echo "No _combined_R2.fastq.gz file found in $dir"
        fi
    fi
done
