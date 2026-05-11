#!/bin/bash

# --- Environment Setup ---
set +e
set -u
ulimit -u 65535 2>/dev/null || true
ulimit -n 65535 2>/dev/null || true

source /mnt/data/personal/sebastian/mambaforge/etc/profile.d/conda.sh

# --- Configuration ---
LOG_DIR="/mnt/data/project/2025_trimmomatic/trimmomatic_paper/logs/scaling_fastp-v1.3.3-no-adapters"
OUT_DIR="/mnt/data/project/2025_trimmomatic/trimmomatic_paper/analysis/scaling_fastp-v1.3.3-no-adapters"
mkdir -p "$LOG_DIR" "$OUT_DIR"


# Input Files
# 1. Plant
PLANT_R1="/mnt/data/project/2025_trimmomatic/trimmomatic_paper/raw/MS01_CasR09_R1_001.fastq.gz"
PLANT_R2="/mnt/data/project/2025_trimmomatic/trimmomatic_paper/raw/MS01_CasR09_R2_001.fastq.gz"

# 2. Human
HUMAN_R1="/mnt/data/project/2025_trimmomatic/trimmomatic_paper/raw/SRR2052337_1.fastq.gz"
HUMAN_R2="/mnt/data/project/2025_trimmomatic/trimmomatic_paper/raw/SRR2052337_2.fastq.gz"

# --- Benchmarking Function ---
run_benchmark() {
    local tool=$1
    local t_count=$2
    local d_name=$3
    local r1_in=$4
    local r2_in=$5

    # Define Output Prefix
    local prefix="${OUT_DIR}/${d_name}_${tool}_${t_count}t"
    local time_log="${LOG_DIR}/${d_name}_${tool}_${t_count}t.time.log"

    echo ">>> Running $tool | Dataset: $d_name | Threads: $t_count"
    
    # Skip if log already exists and is not empty
    if [[ -s "$time_log" ]]; then
        echo "    Skipping (Log exists)."
        return
    fi

    start_time=$(date +%s)

    case $tool in
        "fastp-v1.3.3")
            if conda activate /mnt/bin/fastp/fastp-v1.3.3-conda; then
                /usr/bin/time -v fastp --thread "$t_count" \
                    -i "$r1_in" -I "$r2_in" \
                    -o "${prefix}_R1.fq.gz" -O "${prefix}_R2.fq.gz" \
                    --compression 6 \
                    --json "${prefix}.stress.json" --html "${prefix}.stress.html" \
                    &> "$time_log"
                conda deactivate
            else
                echo "ERROR: fastp env failed."
            fi
            ;;
    esac

    status=$?
    if [[ $status -eq 0 ]]; then
        echo "    Success."
    else
        echo "    FAILED (Exit Code: $status)."
    fi
}

# --- Execution Loops ---

# Threads to test
THREAD_COUNTS=(240 160 80 40 16 8 4 2 1)

# 1. PLANT DATASET LOOP
echo "=========================================================="
echo "STARTING PLANT DATASET"
echo "=========================================================="
for t in "${THREAD_COUNTS[@]}"; do
    run_benchmark "fastp-v1.3.3"       "$t" "Plant" "$PLANT_R1" "$PLANT_R2"
done

# 2. HUMAN DATASET LOOP
echo "=========================================================="
echo "STARTING HUMAN DATASET"
echo "=========================================================="
for t in "${THREAD_COUNTS[@]}"; do
    run_benchmark "fastp-v1.3.3"       "$t" "Human" "$HUMAN_R1" "$HUMAN_R2"
done

echo "All scaling tests complete."