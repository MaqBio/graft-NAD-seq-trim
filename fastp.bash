#!/bin/bash

# Detect if at least 1 sample files were provided
if [ "$#" -eq 0 ]; then
    echo "Usage: $0 sample1 sample2 ..."
    exit 1
fi

# Set up the maximum parallar jobs
MAX_JOBS=3  # 10 threads for each fastp, 30 threads as the maximum
JOB_COUNT=0

# Go through all samplles
for sample in "$@"; do
    # Enter the file
    if [ ! -d "$sample" ]; then
        echo "Warning: Directory $sample does not exist, skipping..."
        continue
    fi
    
    cd "$sample" || exit

    # Set up the file name
    sample_1="${sample}_1.clean.fq.gz"
    sample_2="${sample}_2.clean.fq.gz"
    output_1="${sample}_1.filted.fastq.gz"
    output_2="${sample}_2.filted.fastq.gz"
    report_json="${sample}_report.json"
    report_html="${sample}_report.html"

    # Check the existence of the files
    if [ ! -f "$sample_1" ] || [ ! -f "$sample_2" ]; then
        echo "Warning: Missing files for $sample, skipping..."
        cd ..
        continue
    fi
    
    # Run fastp
    echo "Processing $sample..."
    fastp -i "$sample_1" -I "$sample_2" \
          -l 30 -w 10 -u 20 -n 6 \
          -j "$report_json" -h "$report_html" \
          -o "$output_1" -O "$output_2" &  # & Run the mission parallarly

    ((JOB_COUNT++))

    # If reaching to the maximum threads then wait
    if [ "$JOB_COUNT" -ge "$MAX_JOBS" ]; then
        wait
        JOB_COUNT=0
    fi

    cd ..
done

# Confirm every mission has been done
wait

echo "All selected samples processed."
