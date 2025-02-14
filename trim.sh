#!/bin/bash
set -eu

# set up threasd
THREADS=12

# input the files to be trimmed with space to separate
echo "Please enter file name/.fastq.gz e.g. Sample_SQ24018630-dV1-dV1/SQ24018630-dV1-dV1_combined_R2.fastq.gz. Separate fastq.gz by space"
read -ra FILES_TO_PROCESS

for fq in "${FILES_TO_PROCESS[@]}"; do
    if [[ ! -f "$fq" ]]; then
        echo "Warning: file $fq doesn't exist, skipped"
        continue
    fi

    echo "Processing file: ${fq}"

    # extract reads with barcode GCTTGTTGTG
    fq_branch="${fq%%.fastq.gz}_branch.fastq.gz"
    seqkit grep -j ${THREADS} -s -r -p "^.CTTGTTGT" "$fq" -o "$fq_branch"

    # determine the positions of barcode
    fq_nobarcode="${fq%%.fastq.gz}_nobarcode.fastq.gz"
    seqkit locate -r -p "CTTGTTGT" "$fq_branch" | awk 'NR>1 {print $1"\t"$6+1}' > barcode_positions.txt

    # trim barcode based on their positions
    seqkit subseq -j ${THREADS} --bed barcode_positions.txt "$fq_branch" -o "$fq_nobarcode"

    # remove 3' universal PCR primer
    fq_trimmed="${fq%%.fastq.gz}_3p_trimmed.fastq.gz"
    cutadapt -a "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT" \
        -e 0.2 -m 25 -o "$fq_trimmed" "$fq_nobarcode" 2>&1 | tee -a "${fq}_processing.log"

    echo "Finished processing: ${fq_trimmed}"
done

echo "Finished"
