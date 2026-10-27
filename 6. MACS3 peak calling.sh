#!/usr/bin/env bash
# MACS3 paired-end ATAC-seq peaks; use a filtered BAM and a fresh output directory.
# Supplied genome size: 2,715,809,292 non-mitochondrial A/C/G/T bases (ARS-UCD1.2).
# --keep-dup all retains reads; verify upstream duplicate removal.
macs3 callpeak \
    -t sample.usable.bam -f BAMPE -g 2715809292 -n sample --outdir peaks/sample \
    --keep-dup all -q 0.01 -B --SPMR --call-summits
