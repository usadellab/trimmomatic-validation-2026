# Trimmomatic validation repository

This repository contains the scripts, logs, and analysis code used for the performance evaluation of Trimmomatic v0.40. The study benchmarks Trimmomatic against several state-of-the-art trimming tools.

## Repository structure

```
├── scripts/
│   ├── 01_scaling_benchmark.sh          # Scaling stress test (Plant & Human)
│   ├── 02_residuals_detailed_plant.sh   # Residual adapter assessment (Plant)
│   ├── 03_residuals_detailed_human.sh   # Residual adapter assessment (Human)
│   ├── 04_verify_residuals.py           # Python script for residual verification
│   ├── 05_rabbit_paper_benchmark.sh     # RabbitTrim paper benchmark
│   ├── 06_compression_benchmark.sh      # Parallel vs. Sequential compression test 
│   ├── 07_plot_scaling.R                # R script to visualize scaling
│   └── 08_plot_rabbit_paper_benchmark.R # R script to visualize rabbittrim paper benchmark
│   ├── 09_scaling_benchmark_part2.sh    # Extended scaling stress test (1, 2, 4 threads)
│   ├── 10_plot_scaling_part2.R          # R script to visualize extended scaling
│   ├── 11_filesizes.sh                  # Evaluates compressed output filesizes
│   ├── 12_scaling_benchmark_uncompressed.sh # Scaling benchmark with uncompressed outputs
│   ├── 13_wall_time_uncompressed.sh     # Utility to extract uncompressed wall times
│   ├── 14_scaling_benchmark_fastp-v1.1.0-no-adapters.sh # fastp v1.1.0 scaling (no adapters)
│   ├── 15_scaling_benchmark_fastp-v1.3.3-no-adapters.sh # fastp v1.3.3 scaling (no adapters)
│   ├── 16_scaling_benchmark_fastp-v1.3.3.sh             # fastp v1.3.3 scaling (with adapters)
│   └── 17_plot_comparison_fastp.R       # R script to visualize fastp comparison
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

Three distinct datasets were used to evaluate performance across different species and complexity levels:

1.  **Plant dataset (*Manihot esculenta*)**
    *   **Accession:** ERR15075373
    *   **Description:** Paired-end sequencing data from Cassava.
    *   **Filename in scripts:** `MS01_CasR09_R1_001.fastq.gz`, `MS01_CasR09_R2_001.fastq.gz`

2.  **Human dataset (*Homo sapiens*)**
    *   **Accession:** SRR2052337
    *   **Description:** Paired-end sequencing data from Human.
    *   **Filename in scripts:** `SRR2052337_1.fastq.gz`, `SRR2052337_2.fastq.gz`

3.  **Human dataset (RabbitTrim Benchmark)**
    *   **Accession:** SRR7890824
    *   **Description:** Human dataset used for direct comparison with RabbitTrim publication results.
    *   **Filename in scripts:** `SRR7890824_1.fastq.gz`, `SRR7890824_2.fastq.gz`

## Methodology

The evaluation pipeline consists of the following tests:

**Note**: The shell scripts (`.sh`) are designed for submission to an HPC cluster running Sun Grid Engine (SGE).

### 1. Scaling & Performance Benchmark
**Script:** `01_scaling_benchmark.sh`
*   **Objective:** Comprehensive evaluation of runtime scaling and resource usage across varying thread counts.
*   **Tools:** Trimmomatic, RabbitTrim, fastp, BBDuk, Skewer, Cutadapt.
*   **Datasets:** Runs on both Plant and Human datasets.
*   **Metrics:** Wall-clock time, Peak Memory (RSS), and CPU efficiency.
*   **Configuration:** Tests thread counts of 8, 16, 20, 40, 80, 160, and 240 to determine parallel efficiency limits.

### 2. Accuracy & Residual Assessment
**Scripts:** `02_residuals_detailed_plant.sh`, `03_residuals_detailed_human.sh`, `04_verify_residuals.py`
*   **Objective:** Quantify adapter removal accuracy and detect false positives.
*   **Detection:** Uses `ripgrep` to search for the specific TruSeq3 adapter seed (`AGATCGGAAGAGC`) in the trimmed output files.
*   **Verification:** `04_verify_residuals.py` analyzes the context of residual seeds to distinguish between true adapter remnants and random genomic matches or artifacts.

### 3. RabbitTrim Paper Replication
**Script:** `05_rabbit_paper_benchmark.sh`
*   **Objective:** Replicate the specific benchmark scenario presented in the RabbitTrim publication. Adding the important `compressionLevel` setting to produce comparable results.
*   **Dataset:** SRR7890824 (Human).
*   **Comparison:** Direct comparison of Trimmomatic and RabbitTrim scaling from 8 to 128 threads.

### 4. Compression Efficiency
**Script:** `06_compression_benchmark.sh`
*   **Objective:** Evaluate the performance impact of Trimmomatic's parallel output compression.
*   **Method:** Compares 40-thread parallel compression against single-threaded sequential compression to calculate size overhead and time savings.

### 5. Data Visualization
**Script:** `07_plot_scaling.R`, `08_plot_rabbit_paper_benchmark.R`
*   **Objective:** Visualize the results from the scaling benchmark and residual adapter analysis, aswell as the visualization of the results from the rabbittrim paper benchmark.
*   **Output:** Generates plots for Wall Clock Time, Peak Memory usage, and Fold Speedup across all thread counts and tools. Also displays with the verified residual adapter counts for all tools. Generates Wall Clock Time and Peak Memory usage from the comparison of Trimmomatic with RabbitTrim.

### 6. Extended Evaluation & Uncompressed Benchmarks
**Scripts:** `09_scaling_benchmark_part2.sh`, `10_plot_scaling_part2.R`, `11_filesizes.sh`, `12_scaling_benchmark_uncompressed.sh`, `13_wall_time_uncompressed.sh`
*   **Objective:** Evaluate performance at lower thread counts (1, 2, and 4) and compare the processing speed without the overhead of gzip compression on outputs.
*   **Metrics:** Compression efficiency is also evaluated by comparing the final size of the trimmed reads (`11_filesizes.sh`). Uncompressed runtimes are extracted via the `13_wall_time_uncompressed.sh` utility script.

### 7. Deep fastp evaluation
**Scripts:** `14_scaling_benchmark_fastp-v1.1.0-no-adapters.sh`, `15_scaling_benchmark_fastp-v1.3.3-no-adapters.sh`, `16_scaling_benchmark_fastp-v1.3.3.sh`, `17_plot_comparison_fastp.R` 
* **Objective:** Detailed scaling and accuracy analysis of fastp across different versions (v1.1.0 and v1.3.3) and adapter settings. 
* **Note:** The datapoint for threads = 20 was excluded from this evaluation as fastp v1.3.3 did not terminate when running under this configuration, which is probably due to a software bug. 

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
qsub 01_scaling_benchmark.sh

# 2. Assess accuracy (counts residual adapters)
qsub 02_residuals_detailed_plant.sh
qsub 03_residuals_detailed_human.sh
python3 04_verify_residuals.py analysis/residuals_detailed_plant/* > analysis/residuals_detailed_plant/verified_residuals.txt
python3 04_verify_residuals.py analysis/residuals_detailed_human/* > analysis/residuals_detailed_human/verified_residuals.txt

# 3. Run supplementary rabbittrim paper benchmark
qsub 05_rabbit_paper_benchmark.sh

# 4. Run supplementary compression test
qsub 06_compression_benchmark.sh

# 5. Generate figures
Rscript 07_plot_scaling.R
Rscript 08_plot_rabbit_paper_benchmark.R

# 6. Run extended scaling evaluation and plot
qsub 09_scaling_benchmark_part2.sh
Rscript 10_plot_scaling_part2.R

# 7. Evaluate file sizes and run uncompressed benchmarks
bash 11_filesizes.sh
qsub 12_scaling_benchmark_uncompressed.sh
bash 13_wall_time_uncompressed.sh

# 8. Deep fastp evaluation 
qsub 14_scaling_benchmark_fastp-v1.1.0-no-adapters.sh 
qsub 15_scaling_benchmark_fastp-v1.3.3-no-adapters.sh 
qsub 16_scaling_benchmark_fastp-v1.3.3.sh 

Rscript 17_plot_comparison_fastp.R
```

## License
These scripts are provided under the MIT License to facilitate transparency and reproducibility of the Trimmomatic v0.40 manuscript results.
