#!/bin/bash
#$ -N rabbit_paper_benchmark
#$ -wd /mnt/data/project/2025_trimmomatic/trimmomatic_paper/analysis/
#$ -j y
#$ -o /mnt/data/project/2025_trimmomatic/trimmomatic_paper/logs/rabbit_paper_benchmark.log
#$ -P BigMem.p
#$ -pe smp 128
#$ -l h_vmem=2G
#$ -l h_rt=24:00:00

# --- Environment Setup ---
set +e
set -u
ulimit -u 65535 2>/dev/null || true
ulimit -n 65535 2>/dev/null || true

# --- Configuration ---
LOG_DIR="/mnt/data/project/2025_trimmomatic/trimmomatic_paper/logs/rabbit_scaling"
OUT_DIR="/mnt/data/project/2025_trimmomatic/trimmomatic_paper/analysis/rabbit_scaling"
mkdir -p "$LOG_DIR" "$OUT_DIR"

# Tool Paths
TRIM_JAR="/mnt/bin/trimmomatic/Trimmomatic-0.40/trimmomatic-0.40.jar"
RABBIT_BIN="/mnt/bin/rabbittrim/RabbitTrim_v2.0.0/build/RabbitTrim"
ADAPTERS="/mnt/bin/trimmomatic/Trimmomatic-0.40/adapters/TruSeq3-PE-2-GGGGG.fa"

# Input Files (SRR7890824 - RabbitTrim Paper Dataset)
READ1="/mnt/data/project/2025_trimmomatic/trimmomatic_paper/raw/SRR7890824_1.fastq.gz"
READ2="/mnt/data/project/2025_trimmomatic/trimmomatic_paper/raw/SRR7890824_2.fastq.gz"

# Trimming Steps (Standard Trimmomatic/RabbitTrim equivalent)
TRIM_STEPS="ILLUMINACLIP:${ADAPTERS}:2:30:10 SLIDINGWINDOW:4:20 MINLEN:36"

# --- Benchmarking Function ---
run_benchmark() {
    local tool=$1
    local t_count=$2
    
    # Define Output Prefix
    local prefix="${OUT_DIR}/SRR7890824_${tool}_${t_count}t_c6"
    local time_log="${LOG_DIR}/SRR7890824_${tool}_${t_count}t_c6.time.log"

    echo ">>> Running $tool | Threads: $t_count"
    
    # Skip if log already exists and is not empty
    if [[ -s "$time_log" ]]; then
        echo "    Skipping (Log exists)."
        return
    fi

    case $tool in
        "Trimmomatic")
            /usr/bin/time -v java -jar "$TRIM_JAR" PE -threads "$t_count" -compressLevel 6 \
                "$READ1" "$READ2" \
                "${prefix}_R1.fq.gz" "${prefix}_R1_U.fq.gz" \
                "${prefix}_R2.fq.gz" "${prefix}_R2_U.fq.gz" \
                $TRIM_STEPS \
                &> "$time_log"
            ;;

        "RabbitTrim")
            /usr/bin/time -v "$RABBIT_BIN" trimmomatic \
                --PE \
                --threads "$t_count" \
                --forward "$READ1" \
                --reverse "$READ2" \
                --output "${prefix}_out.gz" \
                --stats "${LOG_DIR}/SRR7890824_${tool}_${t_count}t_c6.stats" \
                --steps $TRIM_STEPS \
                --compressLevel 6 \
                &> "$time_log"       
            ;;
    esac

    status=$?
    if [[ $status -eq 0 ]]; then
        echo "    Success."
    else
        echo "    FAILED (Exit Code: $status)."
    fi
}

# --- Execution Loop ---

# Threads to test (High to Low to fail fast if high count crashes)
THREAD_COUNTS=(128 64 32 16 8)

echo "=========================================================="
echo "STARTING BENCHMARK: SRR7890824 (Human - RabbitTrim Paper)"
echo "Date: $(date)"
echo "=========================================================="

for t in "${THREAD_COUNTS[@]}"; do
    run_benchmark "RabbitTrim"  "$t"
    run_benchmark "Trimmomatic" "$t"
done

echo "Benchmark complete."