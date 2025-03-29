#!/bin/bash

# Ensure there are at least 4 parameters
if [ "$#" -lt 4 ]; then
    echo "Usage: ./readcounts.bash gtf_filepath feature output_filename Sample1 Sample2 ..."
    exit 1
fi

# Analysis
GTF_FILE=$1
FEATURE=$2
OUTPUT_FILE=$3
shift 3  # 移除前 3 个参数，剩下的就是样本名

# Initialize
BAM_FILES=""

# Go through the BAM files in the root file
for SAMPLE in "$@"; do
    SAMPLE_DIR="$SAMPLE" 
    BAM_FILE=$(find "$SAMPLE_DIR" -type f -name "*_sorted.bam")

    if [ -z "$BAM_FILE" ]; then
        echo "Error: No BAM files found for $SAMPLE."
        echo "Please check if the directory exists and contains a *_sorted.bam file."
        exit 1
    else
        BAM_FILES="$BAM_FILES $BAM_FILE"
    fi
done

# Run featureCounts
featureCounts -T 20 -a "$GTF_FILE" -o "$OUTPUT_FILE" -g "$FEATURE" -t exon -p $BAM_FILES
