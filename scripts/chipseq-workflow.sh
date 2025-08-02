#!/bin/bash


set -euo pipefail 

# ChIP-seq data analysis pipeline.
conda activate chip-seq

REF="/home/malm/Tommy/ChIP-seq/data/index/hg19"
Input_control="SRR949142" #SRA ID
Threads=10

# Paths 
DATA_DIR="../data"
RESULTS_DIR="../results"
SCRIPTS_DIR="."

# Create subdirectories
mkdir -p "${RESULTS_DIR}/bams"
mkdir -p "${RESULTS_DIR}/fastqc"
mkdir -p "${RESULTS_DIR}/motifs"
mkdir -p "${RESULTS_DIR}/peaks"

# Check if sra ids textfile exists 
if [[ ! -f "${DATA_DIR}/sra_ids.txt" ]]; then
    echo "sra_ids.txt is not found"
    exit 1
fi

# Process each SRA ID
while read -r SRA; do
    # Skip empty lines and comments
    if [[ -z "$SRA" ]] || [[ "$SRA" == \#* ]]; then
        continue
    fi

    echo "Processing $SRA..."

    # 1 Extract fastq
    fasterq-dump --outdir "${DATA_DIR}" "$SRA"

    # 2 Quality control 
    fastqc -o "${RESULTS_DIR}/fastqc" "${DATA_DIR}/${SRA}.fastq"

    # 3 align with bowtie 
		bowtie2 -x "$REF" -U "${DATA_DIR}/${SRA}.fastq" \
    -S "${DATA_DIR}/${SRA}.sam" \
    -p $Threads --no-unal		

    # 4 SAM to BAM
    samtools view -bS "${DATA_DIR}/${SRA}.sam" | \
    samtools sort -@ $Threads -T "${DATA_DIR}/${SRA}" \
        -o "${RESULTS_DIR}/bams/${SRA}.sorted.bam"

    # 5 index bam
    samtools index "${RESULTS_DIR}/bams/${SRA}.sorted.bam"

    # 6 clean up
    rm "${DATA_DIR}/${SRA}.sam" 
    
    # Only call peaks for non-control samples
    if [[ "$SRA" != "$Input_control" ]]; then
        echo "Calling peaks for $SRA using $Input_control as control..."
        macs2 callpeak \
            -t "${RESULTS_DIR}/bams/${SRA}.sorted.bam" \
            -c "${RESULTS_DIR}/bams/${Input_control}.sorted.bam" \
            -n "${SRA}_peaks" \
            -g hs \
            -p 1e-5 \
            --outdir "${RESULTS_DIR}/peaks"
    fi
done < "${DATA_DIR}/sra_ids.txt"

echo "Workflow completed successfully!"
