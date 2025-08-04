# ChIP-seq Data Processing Pipeline

This is a ChIP-seq processing pipeline, written in bash, that includes the retrieval, quality control, alignment, and peak calling of ChIP-seq datasets using public SRA IDs. The analysis was adopted by the ChIP-seq section in the Biostar Handbook written by Ming Tang. 

https://divingintogeneticsandgenomics.com/files/book_chapters/biostar_handbook_chapter.pdfChIP-seq 

## Overview

The pipeline performs the following steps:

1. **Download** raw FASTQ files using `fasterq-dump`
2. **Quality control** using `FastQC`
3. **Alignment** to the human genome (hg19) using `Bowtie2`
4. **SAM to sorted BAM** conversion with `samtools`
5. **BAM indexing**
6. **Peak calling** using `MACS2` (ChIP vs control)

All steps are performed per sample, using a control sample (specified by `Input_control`) for differential peak calling.

## Requirements

- Bash
- Conda environment with the following tools:
  - `fasterq-dump`
  - `fastqc`
  - `bowtie2`
  - `samtools`
  - `macs2`
- Indexed Bowtie2 reference genome (`hg19`)
- SRA IDs listed in `data/sra_ids.txt`

## Usage

Activate the environment and run the pipeline:

```bash
conda activate chip-seq
./chipseq-workflow.sh

Command being timed: "./chipseq-workflow.sh"
        User time (seconds): 6419.52
        System time (seconds): 210.78
        Percent of CPU this job got: 387%
        Elapsed (wall clock) time (h:mm:ss or m:ss): 28:30.85
        Average shared text size (kbytes): 0
        Average unshared data size (kbytes): 0
        Average stack size (kbytes): 0
        Average total size (kbytes): 0
        Maximum resident set size (kbytes): 7619176
        Average resident set size (kbytes): 0
        Major (requiring I/O) page faults: 1
        Minor (reclaiming a frame) page faults: 3737185
        Voluntary context switches: 6271241
        Involuntary context switches: 8252
        Swaps: 0
        File system inputs: 0
        File system outputs: 8775584
        Socket messages sent: 0
        Socket messages received: 0
        Signals delivered: 0
        Page size (bytes): 4096
        Exit status: 0
```

