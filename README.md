# Trimmomatic validation repository

This repository contains the scripts, logs, and analysis code used for the performance evaluation of Trimmomatic v0.40. The study benchmarks Trimmomatic against several state-of-the-art trimming tools.

## Repository structure

```
├── scripts/
│   ├── 01_tools_eval_plant.sh           # Main benchmark on Manihot esculenta
│   ├── 02_tools_eval_human.sh           # Main benchmark on Homo sapiens
│   ├── 03_plant_accuracy_assessment.sh  # Residual adapter quantification (Plant)
│   ├── 04_human_accuracy_assessment.sh  # Residual adapter quantification (Human)
│   ├── 05_plot_benchmark.R              # R script to generate Figure 3 (Time/Mem/Accuracy)
│   ├── 06_compression_benchmark.sh      # Parallel vs. Sequential compression test
│   ├── 07_scaling_benchmark.sh          # Stress test (20 vs 80 threads)
│   └── 08_plot_scaling_benchmark.R      # R script to visualize scaling efficiency
├── logs/                                # Raw timing and output logs
└── figures/                             # Generated plots (PDF/PNG)
```

## Tools investigated

The following trimming tools were evaluated for performance (runtime/memory) and accuracy (adapter removal):

| Tool | Version | Language | Source |
| :--- | :--- | :--- | :--- |
| **Trimmomatic** | v0.40 | Java | GitHub |
| **RabbitTrim** | v2.0.0 | C++ | Source Build (ISA-L/Zlib linked) |
| **fastp** | v1.1.0 | C++ | Conda |
| **BBDuk** | v39.52 | Java | Conda |
| **Skewer** | v0.2.2 | C++ | Conda |
| **Cutadapt** | v5.2 | Python/C | Conda |

## Datasets

Two distinct datasets were used to evaluate performance across different species and complexity levels:

1.  **Plant dataset (*Manihot esculenta*)**
    *   **Accession:** ERR15075373
    *   **Description:** Paired-end sequencing data from Cassava.
    *   **Filename in scripts:** `MS01_CasR09_R1_001.fastq.gz`, `MS01_CasR09_R2_001.fastq.gz`

2.  **Human dataset (*Homo sapiens*)**
    *   **Accession:** SRR2052337
    *   **Description:** Paired-end sequencing data from Human.
    *   **Filename in scripts:** `SRR2052337_1.fastq.gz`, `SRR2052337_2.fastq.gz`

## Methodology

The evaluation pipeline consists of the following tests:

**Note**: The shell scripts (`.sh`) are designed for submission to an HPC cluster running Sun Grid Engine (SGE).

### 1. Performance benchmarking
**Scripts:** `01_tools_eval_plant.sh`, `02_tools_eval_human.sh`
*   Measures wall-clock time, CPU usage, and peak memory using `/usr/bin/time -v`.
*   Runs all tools with a standardized thread count (40 threads) on the same HPC cluster.
*   Standardized trimming parameters are used where possible (e.g., adapter clipping, sliding window quality trimming).
*   Parameters: Standard TruSeq3 adapter trimming with a sliding window quality filter (Window:4, Quality:20) and minimum length of 36bp.

### 2. Accuracy assessment
**Scripts:** `03_plant_accuracy_assessment.sh`, `04_human_accuracy_assessment.sh`
*   Evaluates the effectiveness of adapter removal.
*   Searches for residual adapter seeds (`AGATCGGAAGAGC`) using parallel decompression (`pragzip`) piped to SIMD-accelerated grep (`ripgrep`).

### 3. Scaling stress test
**Script:** `07_scaling_benchmark.sh`
*   Evaluates how well the tools scale with available resources.
*   Runs benchmarks at varying thread counts (20 vs 80 threads) to observe speedup efficiency.

### 4. Compression efficiency
**Script:** `06_compression_benchmark.sh`
*   Specifically tests the parallel output compression feature in Trimmomatic v0.40.
*   Compares file size and compression time between parallel (multi-threaded) and sequential (single-threaded) gzip output modes.

### 5. Data visualization
**Scripts:** `05_plot_benchmark.R`, `08_plot_scaling_benchmark.R`
*   **Benchmark Summary:** `05_plot_benchmark.R` parses the timing logs and accuracy counts to generate bar charts comparing Wall Clock Time, Peak Memory, and Residual Adapters.
*   **Scaling Analysis:** `08_plot_scaling_benchmark.R` visualizes the scaling stress test results, plotting runtime against thread count to demonstrate parallel efficiency.
*   **Output:** All generated plots are saved to the `figures/` directory.

## Validation environment

| Component | Details | 
| :--- | :--- | 
| OS | `Ubuntu 22.04.5 LTS` |
| Kernel | `5.15.0-143-generic` |
| Model | `AMD EPYC 7763 64-Core Processor` |
| CPUs | `256` |
| Sockets | `2` |
| Threads per core | `2` |
| Cores per socket | `64` |
| RAM | `2 Tb` |

## Reproducing results
The scripts are numbered in the order they should be executed.
```
# 1. Run main benchmarks (ensure input paths in scripts match your system)
qsub 01_tools_eval_plant.sh
qsub 02_tools_eval_human.sh

# 2. Assess accuracy (counts residual adapters)
qsub 03_plant_accuracy_assessment.sh
qsub 04_human_accuracy_assessment.sh

# 3. Generate main comparison figure
Rscript 05_plot_benchmark.R

# 4. Run supplementary compression test
qsub 06_compression_benchmark.sh

# 5. Run scaling stress test
qsub 07_scaling_benchmark.sh

# 6. Generate scaling figure
Rscript 08_plot_scaling_benchmark.R
```

## License
These scripts are provided under the MIT License to facilitate transparency and reproducibility of the Trimmomatic v0.40 manuscript results.
