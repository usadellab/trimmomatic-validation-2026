#!/bin/bash
#$ -N scaling_stress_test
#$ -wd /mnt/data/project/2025_trimmomatic/trimmomatic_paper/analysis/
#$ -j y
#$ -o /mnt/data/project/2025_trimmomatic/trimmomatic_paper/logs/scaling_stress_test.log
#$ -pe smp 80
#$ -l h_vmem=2G
#$ -l h_rt=24:00:00

# --- Environment Setup ---
set +e
set -u
ulimit -u 65535 2>/dev/null || true
ulimit -n 65535 2>/dev/null || true

source /mnt/data/personal/mambaforge/etc/profile.d/conda.sh

# --- Configuration ---
LOG_DIR="/mnt/data/project/2025_trimmomatic/trimmomatic_paper/logs/scaling"
OUT_DIR="/mnt/data/project/2025_trimmomatic/trimmomatic_paper/analysis/scaling"
mkdir -p "$LOG_DIR" "$OUT_DIR"

# Tool Paths
TRIM_JAR="/mnt/bin/trimmomatic/Trimmomatic-0.40/trimmomatic-0.40.jar"
RABBIT_BIN="/mnt/bin/rabbittrim/RabbitTrim_v2.0.0/build/RabbitTrim"
ADAPTERS="/mnt/bin/trimmomatic/Trimmomatic-0.40/adapters/TruSeq3-PE-2-GGGGG.fa"
TRIM_STEPS="ILLUMINACLIP:${ADAPTERS}:2:30:10 SLIDINGWINDOW:4:20 MINLEN:36"

# Input Files
# 1. Plant data
PLANT_R1="/mnt/data/project/2025_trimmomatic/trimmomatic_paper/raw/MS01_CasR09_R1_001.fastq.gz"
PLANT_R2="/mnt/data/project/2025_trimmomatic/trimmomatic_paper/raw/MS01_CasR09_R2_001.fastq.gz"

# 2. Human data
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
        "Trimmomatic")
            /usr/bin/time -v java -jar "$TRIM_JAR" PE -threads "$t_count" \
                "$r1_in" "$r2_in" \
                "${prefix}_R1.fq.gz" "${prefix}_R1_U.fq.gz" \
                "${prefix}_R2.fq.gz" "${prefix}_R2_U.fq.gz" \
                ILLUMINACLIP:"$ADAPTERS":2:30:10 SLIDINGWINDOW:4:20 MINLEN:36 \
                &> "$time_log"
            ;;

        "RabbitTrim")
            # RabbitTrim usually needs specific flags, using Trimmomatic mode
            /usr/bin/time -v "$RABBIT_BIN" trimmomatic \
                --PE \
                --threads "$t_count" \
                --forward "$r1_in" \
                --reverse "$r2_in" \
                --output "${prefix}_out" \
                --stats "${LOG_DIR}/${d_name}_${tool}_${t_count}t.stats" \
                --steps $TRIM_STEPS \
                --compressLevel 1 \
                &> "$time_log"
            ;;

        "fastp")
            if conda activate /mnt/bin/fastp/fastp-v1.1.0-conda; then
                /usr/bin/time -v fastp --thread "$t_count" \
                    -i "$r1_in" -I "$r2_in" \
                    -o "${prefix}_R1.fq.gz" -O "${prefix}_R2.fq.gz" \
                    --adapter_fasta "$ADAPTERS" \
                    --json "${prefix}.stress.json" --html "${prefix}.stress.html" \
                    &> "$time_log"
                conda deactivate
            else
                echo "ERROR: fastp env failed."
            fi
            ;;

        "BBDuk")
            if conda activate /mnt/bin/bbmap/bbmap-v39.52-conda; then
                # Adjust memory for BBDuk (safe limit 20G to avoid SGE kill)
                /usr/bin/time -v bbduk.sh -Xmx20g threads="$t_count" ordered=f \
                    in="$r1_in" in2="$r2_in" \
                    out="${prefix}_R1.fq.gz" out2="${prefix}_R2.fq.gz" \
                    ref="$ADAPTERS" ktrim=r k=23 mink=11 hdist=1 tpe tbo \
                    &> "$time_log"
                conda deactivate
            else
                echo "ERROR: BBDuk env failed."
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
THREAD_COUNTS=(20 80)

# 1. PLANT DATASET LOOP
echo "=========================================================="
echo "STARTING PLANT DATASET"
echo "=========================================================="
for t in "${THREAD_COUNTS[@]}"; do
    run_benchmark "Trimmomatic" "$t" "Plant" "$PLANT_R1" "$PLANT_R2"
    run_benchmark "fastp"       "$t" "Plant" "$PLANT_R1" "$PLANT_R2"
    run_benchmark "RabbitTrim"  "$t" "Plant" "$PLANT_R1" "$PLANT_R2"
    run_benchmark "BBDuk"       "$t" "Plant" "$PLANT_R1" "$PLANT_R2"
done

# 2. HUMAN DATASET LOOP
echo "=========================================================="
echo "STARTING HUMAN DATASET"
echo "=========================================================="
for t in "${THREAD_COUNTS[@]}"; do
    run_benchmark "Trimmomatic" "$t" "Human" "$HUMAN_R1" "$HUMAN_R2"
    run_benchmark "fastp"       "$t" "Human" "$HUMAN_R1" "$HUMAN_R2"
    run_benchmark "RabbitTrim"  "$t" "Human" "$HUMAN_R1" "$HUMAN_R2"
    run_benchmark "BBDuk"       "$t" "Human" "$HUMAN_R1" "$HUMAN_R2"
done

echo "All scaling tests complete."