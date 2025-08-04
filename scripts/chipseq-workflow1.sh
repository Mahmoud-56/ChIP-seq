#!/bin/bash
set -euo pipefail

conda activate chip-seq

REF="/home/malm/Tommy/ChIP-seq/data/index/hg19"
Input_control="SRR949142"
Threads=10

# Paths
DATA_DIR="../data"
RESULTS_DIR="../results"
SCRIPTS_DIR="."

# Create subdirectories
mkdir -p "${RESULTS_DIR}/bams" "${RESULTS_DIR}/fastqc" "${RESULTS_DIR}/motifs" "${RESULTS_DIR}/peaks" "${RESULTS_DIR}/logs"

# Check if sra ids textfile exists
if [[ ! -f "${DATA_DIR}/sra_ids.txt" ]]; then
    echo "sra_ids.txt is not found"
    exit 1
fi

# Process control sample first
echo "Processing control sample $Input_control..."

fasterq-dump --outdir "${DATA_DIR}" "$Input_control"

fastqc -o "${RESULTS_DIR}/fastqc" "${DATA_DIR}/${Input_control}.fastq"

bowtie2 -x "$REF" -U "${DATA_DIR}/${Input_control}.fastq" \
    -S "${DATA_DIR}/${Input_control}.sam" \
    -p "$Threads" --no-unal

samtools view -bS "${DATA_DIR}/${Input_control}.sam" | \
    samtools sort -@ "$Threads" -T "${DATA_DIR}/${Input_control}" \
    -o "${RESULTS_DIR}/bams/${Input_control}.sorted.bam"

samtools index "${RESULTS_DIR}/bams/${Input_control}.sorted.bam"

rm "${DATA_DIR}/${Input_control}.sam"

echo "Parallel processing IP samples..."
tail -n +2 "${DATA_DIR}/sra_ids.txt" | awk '{$1=$1; print}' | parallel -j "$Threads" --joblog "${RESULTS_DIR}/logs/parallel.log" '
    echo "Processing {}..."

    # 1. FASTQ extraction
    fasterq-dump --outdir '"${DATA_DIR}"' {}

    # 2. Quality control
    fastqc -o '"${RESULTS_DIR}"'/fastqc '"${DATA_DIR}"'/{}.fastq

    # 3. Alignment
    bowtie2 -x '"${REF}"' -U '"${DATA_DIR}"'/{}.fastq \
        -S '"${DATA_DIR}"'/{}.sam \
        -p 1 --no-unal

    # 4. SAM to BAM
    samtools view -bS '"${DATA_DIR}"'/{}.sam | \
        samtools sort -@ 1 -T '"${DATA_DIR}"'/{} \
        -o '"${RESULTS_DIR}"'/bams/{}.sorted.bam

    # 5. Index BAM
    samtools index '"${RESULTS_DIR}"'/bams/{}.sorted.bam

    # 6. Cleanup
    rm '"${DATA_DIR}"'/{}.sam

    # 7. Peak calling
    if [[ -f '"${RESULTS_DIR}"'/bams/'"${Input_control}"'.sorted.bam ]]; then
        macs2 callpeak \
            -t '"${RESULTS_DIR}"'/bams/{}.sorted.bam \
            -c '"${RESULTS_DIR}"'/bams/'"${Input_control}"'.sorted.bam \
            -n {}_peaks \
            -g hs \
            -p 1e-5 \
            --outdir '"${RESULTS_DIR}"'/peaks
    else
        echo "Control BAM missing for {}" >&2
        exit 1
    fi
' &> "${RESULTS_DIR}/logs/parallel_processing.log"

