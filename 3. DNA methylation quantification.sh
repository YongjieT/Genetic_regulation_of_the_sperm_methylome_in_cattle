#!/usr/bin/env bash
# Staged WGBS commands; configure working directories and SAMPLE.
set -euo pipefail
: "${SAMPLE:?Set SAMPLE}"

# Quality control
cd ~/QC
trim_galore --paired --fastqc -q 25 --length 50 -o ~/QC \
${SAMPLE}/${SAMPLE}_1.fq.gz ${SAMPLE}/${SAMPLE}_2.fq.gz

# Alignment
cd ~/Align
bismark --bowtie2 -o ~/Align \
--multicore 4 -p 4 -N 1 ~/1REFpreparation/ \
-1 ~/QC/${SAMPLE}_1_val_1.fq.gz \
-2 ~/QC/${SAMPLE}_2_val_2.fq.gz

# Deduplication
cd ~/deduplicate_bismark
deduplicate_bismark -p --bam ~/Align/${SAMPLE}_1_val_1_bismark_bt2_pe.bam

# Initialize the reference once; supply sorted/indexed BAMs before bam2pat.
wgbstools init_genome ARSUCD1.2 --fasta_path Bos_taurus.ARS-UCD1.2.dna.toplevel.fa
wgbstools bam2pat ${SAMPLE}_1_val_1_bismark_bt2_pe.deduplicated_sorted.bam
# Segment the intended cohort BETA files; require at least two CpGs per block.
wgbstools segment --betas *.beta --min_cpg 2 --max_bp 5000 -o 125_cattle_blocks.small.bed
wgbstools index 125_cattle_blocks.small.bed
wgbstools beta_to_table 125_cattle_blocks.small.bed -c 10 --betas *.beta > 125_cattle_DNAm_blocks_level.txt
