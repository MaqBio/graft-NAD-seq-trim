#!/bin/bash

# Check the parameters
if [ "$#" -lt 3 ]; then
    echo "Usage: $0 index gtf sample1 sample2 ..."
    exit 1
fi

# Analyze the parameters
INDEX=$1
GTF=$2
shift 2  # Shift the parameters, make sure the parameters after the first 2 are samples

# Set up 20 threads for each step
THREADS=20

# Determine f(x) for mapping
run_hisat2() {
    SAMPLE=$1
    INDEX=$2

    # Enter the directory
    if [ ! -d "$SAMPLE" ]; then
        echo "Warning: Directory $SAMPLE does not exist, skipping..."
        return 1
    fi
    cd "$SAMPLE" || return 1

    # Define the files
    SAMPLE_1="${SAMPLE}_1.filted.fastq.gz"
    SAMPLE_2="${SAMPLE}_2.filted.fastq.gz"
    SAM_OUTPUT="${SAMPLE}.sam"
    BAM_OUTPUT="${SAMPLE}.bam"
    SORTED_BAM_OUTPUT="${SAMPLE}_sorted.bam"
    LOG_FILE="${SAMPLE}_hisat2map.log"

    # Check if filtered.fastq.gz files exist
    if [ ! -f "$SAMPLE_1" ] || [ ! -f "$SAMPLE_2" ]; then
        echo "Warning: Missing fastq files for $SAMPLE, skipping..."
        cd ..
        return 1
    fi

    # Run HISAT2 alignment
    echo "Running HISAT2 for $SAMPLE..."
    hisat2 -p "$THREADS" -x "$INDEX" -1 "$SAMPLE_1" -2 "$SAMPLE_2" -S "$SAM_OUTPUT" 2> "$LOG_FILE"

    # Convert SAM to BAM
    echo "Converting SAM to BAM for $SAMPLE..."
    samtools view -@ "$THREADS" -bS "$SAM_OUTPUT" > "$BAM_OUTPUT"

    # Sort BAM files
    echo "Sorting BAM file for $SAMPLE..."
    samtools sort -@ "$THREADS" -o "$SORTED_BAM_OUTPUT" "$BAM_OUTPUT"

    # Clean up the original SAM/BAM (optional)
    rm -f "$SAM_OUTPUT" "$BAM_OUTPUT"

    echo "Finished processing $SAMPLE."
    cd ..
}

# Process each sample sequentially
echo "Starting sequential HISAT2 alignment..."
for SAMPLE in "$@"; do
    run_hisat2 "$SAMPLE" "$INDEX"
done

# Check if every file has the corresponding _sorted.bam
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

# If there are failed samples, stop featureCounts
if [ "${#FAILED_SAMPLES[@]}" -gt 0 ]; then
    echo "Some samples failed. Skipping featureCounts to avoid incomplete data."
    exit 1
fi

# Run featureCounts on the sorted BAM files
echo "Running featureCounts on sorted BAM files..."
featureCounts -F GTF -C -T 20 -M -t gene -g gene_id -O -p -a "$GTF" \
    -o merged.txt "${SORTED_BAM_FILES[@]}"

echo "All samples processed successfully."
