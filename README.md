# Trimmomatic Benchmarking Repository

This repository contains the scripts, logs, and analysis code used for the performance evaluation of Trimmomatic v0.40. The study benchmarks Trimmomatic against several state-of-the-art trimming tools.

## Repository Structure

```
├── scripts/
│   ├── 01_tools_eval_plant.sh        # Main benchmark on Manihot esculenta
│   ├── 02_tools_eval_human.sh        # Main benchmark on Homo sapiens
│   ├── 03_plant_accuracy_assessment.sh  # Residual adapter quantification (Plant)
│   ├── 04_human_accuracy_assessment.sh  # Residual adapter quantification (Human)
│   ├── 05_plot_benchmark.R           # R script to generate Figure 3 (Time/Mem/Accuracy)
│   ├── 06_compression_benchmark.sh   # Parallel vs. Sequential compression test
│   ├── 07_scaling_benchmark.sh       # Stress test (20 vs 80 threads)
│   └── 08_plot_scaling_benchmark.R   # R script to visualize scaling efficiency
├── logs/                             # Raw timing and output logs
└── figures/                          # Generated plots (PDF/PNG)
```

## Tools Investigated

The following trimming tools were evaluated for performance (runtime/memory) and accuracy (adapter removal):

1.  **Trimmomatic** (Java)
2.  **fastp** (C++)
3.  **BBDuk** (Java / BBMap suite)
4.  **Skewer** (C++)
5.  **Cutadapt** (Python)
6.  **RabbitTrim** (C++)

## Datasets

Two distinct datasets were used to evaluate performance across different species and complexity levels:

1.  **Plant Dataset (*Manihot esculenta*)**
    *   **Accession:** ERR15075373
    *   **Description:** Paired-end sequencing data from Cassava.
    *   **Filename in scripts:** `MS01_CasR09_R1_001.fastq.gz`, `MS01_CasR09_R2_001.fastq.gz`

2.  **Human Dataset (*Homo sapiens*)**
    *   **Accession:** SRR2052337
    *   **Description:** Paired-end sequencing data from Human.
    *   **Filename in scripts:** `SRR2052337_1.fastq.gz`, `SRR2052337_2.fastq.gz`

## Methodology

The evaluation pipeline consists of the following tests:

### 1. Performance Benchmarking
**Scripts:** `01_tools_eval_plant.sh`, `02_tools_eval_human.sh`
*   Measures wall-clock time, CPU usage, and peak memory using `/usr/bin/time -v`.
*   Runs all tools with a standardized thread count (40 threads) on the same HPC cluster.
*   Standardized trimming parameters are used where possible (e.g., adapter clipping, sliding window quality trimming).

### 2. Accuracy Assessment
**Scripts:** `03_plant_accuracy_assessment.sh`, `04_human_accuracy_assessment.sh`
*   Evaluates the effectiveness of adapter removal.
*   Searches for residual adapter seeds (e.g., `AGATCGGAAGAGC`) in the trimmed output files using `pragzip` (`rapidgzip`) and `grep` or `ripgrep`.

### 3. Scaling Stress Test
**Script:** `07_scaling_benchmark.sh`
*   Evaluates how well the tools scale with available resources.
*   Runs benchmarks at varying thread counts (e.g., 20 vs 80 threads) to observe speedup efficiency.

### 4. Compression Efficiency
**Script:** `06_compression_benchmark.sh`
*   Specifically tests the parallel output compression feature in Trimmomatic v0.40.
*   Compares file size and compression time between parallel (multi-threaded) and sequential (single-threaded) gzip output modes.

### 5. Data Visualization
**Scripts:** `05_plot_benchmark.R`, `08_plot_scaling_benchmark.R`
*   **Benchmark Summary:** `05_plot_benchmark.R` parses the timing logs and accuracy counts to generate bar charts comparing Wall Clock Time, Peak Memory, and Residual Adapters (Figure 3).
*   **Scaling Analysis:** `08_plot_scaling_benchmark.R` visualizes the scaling stress test results, plotting runtime against thread count to demonstrate parallel efficiency.
*   **Output:** All generated plots are saved to the `figures/` directory.
