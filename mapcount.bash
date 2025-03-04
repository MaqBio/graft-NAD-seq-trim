#!/bin/bash

# Check the parameters
if [ "$#" -lt 3 ]; then
    echo "Usage: $0 index gtf sample1 sample2 ..."
    exit 1
fi

# Analyze the parameters
INDEX=$1
GTF=$2
shift 2  # Shift the parameters, make sure the parameters after the first 2 were samples

# Set up the max threads
THREADS_PER_SAMPLE=10
MAX_JOBS=2  # Two samples run in parallel

# Determine f(x) for mapping
run_hisat2() {
    SAMPLE=$1
    INDEX=$2

    # Enter the file
    if [ ! -d "$SAMPLE" ]; then
        echo "Warning: Directory $SAMPLE does not exist, skipping..."
        return 1
    fi
    cd "$SAMPLE" || return 1

    # Define the file
    SAMPLE_1="${SAMPLE}_1.filted.fastq.gz"
    SAMPLE_2="${SAMPLE}_2.filted.fastq.gz"
    SAM_OUTPUT="${SAMPLE}.sam"
    BAM_OUTPUT="${SAMPLE}.bam"
    SORTED_BAM_OUTPUT="${SAMPLE}_sorted.bam"
    LOG_FILE="${SAMPLE}_hisat2map.log"

    # Check if filtered.fastq.gz file exist
    if [ ! -f "$SAMPLE_1" ] || [ ! -f "$SAMPLE_2" ]; then
        echo "Warning: Missing fastq files for $SAMPLE, skipping..."
        cd ..
        return 1
    fi

    # Hisat2 alignment
    echo "Running HISAT2 for $SAMPLE..."
    hisat2 -p "$THREADS_PER_SAMPLE" -x "$INDEX" -1 "$SAMPLE_1" -2 "$SAMPLE_2" -S "$SAM_OUTPUT" 2> "$LOG_FILE"

    # Convert SAM to BAM
    echo "Converting SAM to BAM for $SAMPLE..."
    samtools view -@ "$THREADS_PER_SAMPLE" -bS "$SAM_OUTPUT" > "$BAM_OUTPUT"

    # Sort BAM files
    echo "Sorting BAM file for $SAMPLE..."
    samtools sort -@ "$THREADS_PER_SAMPLE" -o "$SORTED_BAM_OUTPUT" "$BAM_OUTPUT"

    # Clean up the original SAM/BAM (optional)
    rm -f "$SAM_OUTPUT" "$BAM_OUTPUT"

    echo "Finished processing $SAMPLE."
    cd ..
}

export -f run_hisat2

# Run HISAT2+SAMtools for different samples in parallel. We will limit parallel jobs to 2 at a time.
echo "Starting parallel HISAT2 alignment..."
parallel -j "$MAX_JOBS" run_hisat2 {} "$INDEX" ::: "$@"

# Check if every file has the corresponding `_sorted.bam`
echo "Checking sorted BAM files..."
FAILED_SAMPLES=()
SORTED_BAM_FILES=()
for SAMPLE in "$@"; do
    if [ -f "${SAMPLE}/${SAMPLE}_sorted.bam" ]; then
        SORTED_BAM_FILES+=("${SAMPLE}/${SAMPLE}_sorted.bam")
    else
        echo "Error: ${SAMPLE}_sorted.bam is missing!"
        FAILED_SAMPLES+=("$SAMPLE")
    fi
done

# If there is a lack of sorted.bam then stop featureCounts
if [ "${#FAILED_SAMPLES[@]}" -gt 0 ]; then
    echo "Some samples failed. Skipping featureCounts to avoid incomplete data."
    exit 1
fi

# Run featureCounts
echo "Running featureCounts on sorted BAM files..."
featureCounts -F GTF -C -T "$TOTAL_THREADS" -M -t gene -g gene_id -O -s 1 -a "$GTF" \
    -o merged.txt --fraction "${SORTED_BAM_FILES[@]}"

echo "All samples processed successfully."
