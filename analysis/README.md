# Planarian RNA-seq pipeline for GSE171247 (wtap RNAi)

This folder contains a lightweight, reproducible RNA-seq pipeline tailored to the GSE171247 bulk RNA-seq subseries (planarian regeneration following **wtap** RNAi). The steps prioritize transparent preprocessing, quantification, and differential expression with **DESeq2**. Update the metadata file before running to ensure sample annotations match your downloads.

## Data inputs
- **SRR accessions (initial focus):** SRR14119237, SRR14119220, SRR14119230, SRR14119224, SRR14119234, SRR14119222.
- Place a reference transcriptome FASTA for *Schmidtea mediterranea* at `refs/planarian_transcriptome.fa` and build the Salmon index under `refs/salmon_index/` (see below).
- Edit `metadata_samples.csv` to reflect the correct condition and time point for each SRR before running DE analysis.

## Quickstart
1. Create a conda environment with the required tools:
   ```bash
   mamba create -n planaria-rnaseq fastp salmon sra-tools multiqc -c conda-forge -c bioconda
   mamba activate planaria-rnaseq
   ```
2. Download reads, trim adapters, and quantify with Salmon:
   ```bash
   bash scripts/download_and_quant.sh
   ```
3. Summarize quantifications to gene-level counts and run DESeq2:
   ```bash
   Rscript scripts/deseq2_analysis.R
   ```

## Reference indexing (once per reference)
```bash
mkdir -p refs
wget -O refs/planarian_transcriptome.fa "<URL-to-transcriptome>"
salmon index -t refs/planarian_transcriptome.fa -i refs/salmon_index --type quasi -k 31
```

## Files
- `metadata_samples.csv` — sample annotations for the six SRR accessions (condition, time point, replicate).
- `scripts/download_and_quant.sh` — download, QC (fastp), and Salmon quasi-mapping/quantification.
- `scripts/deseq2_analysis.R` — DESeq2 workflow using `tximport` on Salmon quantifications with visualization stubs (PCA, volcano).

## Outputs
After running the scripts, expect:
- `fastq/` — raw FASTQ files from `fasterq-dump`.
- `qc/fastp/` — HTML/QC reports from fastp.
- `salmon/` — Salmon quantifications (one folder per SRR).
- `results/deseq2/` — DE results tables and plots (PCA, volcano, heatmap).

## Notes
- The pipeline assumes paired-end reads; adjust the flags in `download_and_quant.sh` if single-end libraries are detected.
- Replace the transcript-to-gene mapping section in `deseq2_analysis.R` with a GTF/annotation appropriate for your transcriptome build.
- MultiQC can be added after fastp/Salmon to aggregate QC reports.
