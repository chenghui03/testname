#!/usr/bin/env bash
set -euo pipefail

# Download, QC, and quantify the six SRR runs for GSE171247.
# Requirements: sra-tools (prefetch/fasterq-dump), fastp, salmon.

mkdir -p fastq qc/fastp salmon logs

# Download and convert SRA to FASTQ (paired-end assumed).
for srr in SRR14119237 SRR14119220 SRR14119230 SRR14119224 SRR14119234 SRR14119222; do
  if [ ! -f fastq/${srr}_1.fastq.gz ]; then
    echo "[INFO] Fetching ${srr}" | tee -a logs/download.log
    prefetch --max-size 200G ${srr}
    fasterq-dump --split-files --threads 8 --outdir fastq ${srr}
    gzip fastq/${srr}_1.fastq fastq/${srr}_2.fastq
  else
    echo "[SKIP] FASTQ for ${srr} already present" | tee -a logs/download.log
  fi

done

echo "[INFO] Running fastp for adapter/quality trimming" | tee -a logs/fastp.log
for srr in SRR14119237 SRR14119220 SRR14119230 SRR14119224 SRR14119234 SRR14119222; do
  fastp \
    -i fastq/${srr}_1.fastq.gz -I fastq/${srr}_2.fastq.gz \
    -o fastq/${srr}_1.trim.fastq.gz -O fastq/${srr}_2.trim.fastq.gz \
    --thread 8 \
    --html qc/fastp/${srr}.fastp.html \
    --json qc/fastp/${srr}.fastp.json \
    --detect_adapter_for_pe

done

echo "[INFO] Salmon quantification" | tee -a logs/salmon.log
for srr in SRR14119237 SRR14119220 SRR14119230 SRR14119224 SRR14119234 SRR14119222; do
  salmon quant -i refs/salmon_index -l A \
    -1 fastq/${srr}_1.trim.fastq.gz -2 fastq/${srr}_2.trim.fastq.gz \
    -o salmon/${srr} --threads 8 --validateMappings --gcBias

done
