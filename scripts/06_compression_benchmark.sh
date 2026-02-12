#!/bin/bash
#$ -N trimmomatic_compression_test
#$ -wd /mnt/data/project/2025_trimmomatic/trimmomatic_paper/analysis/
#$ -j y
#$ -o /mnt/data/project/2025_trimmomatic/trimmomatic_paper/logs/compression_eval.log
#$ -pe smp 40
#$ -l h_vmem=1G
#$ -l h_rt=12:00:00

# --- Environment Setup ---
set -u
set -e

# --- Configuration ---
RAW_DIR="/mnt/data/project/2025_trimmomatic/trimmomatic_paper/raw"
OUT_DIR="/mnt/data/project/2025_trimmomatic/trimmomatic_paper/analysis/compression_test"
LOG_DIR="/mnt/data/project/2025_trimmomatic/trimmomatic_paper/logs"
mkdir -p "$OUT_DIR" "$LOG_DIR"

# Input Files (Plant Data)
READ1="${RAW_DIR}/MS01_CasR09_R1_001.fastq.gz"
READ2="${RAW_DIR}/MS01_CasR09_R2_001.fastq.gz"

# Trimmomatic Settings
TRIMMOMATIC_JAR="/mnt/bin/trimmomatic/Trimmomatic-0.40/trimmomatic-0.40.jar"
ADAPTERS="/mnt/bin/trimmomatic/Trimmomatic-0.40/adapters/TruSeq3-PE-2-GGGGG.fa"
TRIM_STEPS="ILLUMINACLIP:${ADAPTERS}:2:30:10 SLIDINGWINDOW:4:20 MINLEN:36"

echo "=== Compression Benchmark: Plant ==="
echo "Host: $(hostname)"
echo "Date: $(date)"

# --- 1. Parallel Compression (40 Threads) ---
echo "----------------------------------------------------------------"
echo "Starting PARALLEL Run (40 Threads)..."
start_time=$(date +%s)

/usr/bin/time -v java -jar "$TRIMMOMATIC_JAR" PE -threads 40 \
    "$READ1" "$READ2" \
    "$OUT_DIR/parallel_R1.fq.gz" "$OUT_DIR/parallel_R1_unpaired.fq.gz" \
    "$OUT_DIR/parallel_R2.fq.gz" "$OUT_DIR/parallel_R2_unpaired.fq.gz" \
    $TRIM_STEPS \
    2> "$LOG_DIR/compression_parallel.log"

end_time=$(date +%s)
echo "Parallel Run finished in $((end_time - start_time)) seconds."

# --- 2. Sequential Compression (1 Thread) ---
echo "----------------------------------------------------------------"
echo "Starting SEQUENTIAL Run (1 Thread)..."
# NOTE: This effectively forces single-threaded compression
start_time=$(date +%s)

/usr/bin/time -v java -jar "$TRIMMOMATIC_JAR" PE -threads 1 \
    "$READ1" "$READ2" \
    "$OUT_DIR/sequential_R1.fq.gz" "$OUT_DIR/sequential_R1_unpaired.fq.gz" \
    "$OUT_DIR/sequential_R2.fq.gz" "$OUT_DIR/sequential_R2_unpaired.fq.gz" \
    $TRIM_STEPS \
    2> "$LOG_DIR/compression_sequential.log"

end_time=$(date +%s)
echo "Sequential Run finished in $((end_time - start_time)) seconds."

# --- 3. Size Comparison Report ---
echo "----------------------------------------------------------------"
echo "COMPRESSION RATIO EVALUATION"
echo "----------------------------------------------------------------"
printf "%-20s | %-15s | %-15s\n" "Mode" "R1 Size (Bytes)" "R2 Size (Bytes)"
echo "---------------------|-----------------|-----------------"

# Get file sizes in bytes
size_par_r1=$(stat -c%s "$OUT_DIR/parallel_R1.fq.gz")
size_par_r2=$(stat -c%s "$OUT_DIR/parallel_R2.fq.gz")
size_seq_r1=$(stat -c%s "$OUT_DIR/sequential_R1.fq.gz")
size_seq_r2=$(stat -c%s "$OUT_DIR/sequential_R2.fq.gz")

printf "%-20s | %-15s | %-15s\n" "Parallel (40t)" "$size_par_r1" "$size_par_r2"
printf "%-20s | %-15s | %-15s\n" "Sequential (1t)" "$size_seq_r1" "$size_seq_r2"

echo "---------------------|-----------------|-----------------"

# Calculate percentage difference
# Formula: ((Parallel - Sequential) / Sequential) * 100
diff_r1=$(awk "BEGIN {print (($size_par_r1 - $size_seq_r1) / $size_seq_r1) * 100}")
diff_r2=$(awk "BEGIN {print (($size_par_r2 - $size_seq_r2) / $size_seq_r2) * 100}")

echo "Size Difference (Parallel vs Sequential):"
echo "R1: $diff_r1 %"
echo "R2: $diff_r2 %"

echo "----------------------------------------------------------------"
echo "Done."