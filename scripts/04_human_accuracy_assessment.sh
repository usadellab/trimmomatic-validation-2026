#!/bin/bash
#$ -N human_adapter_assessment
#$ -wd /mnt/data/project/2025_trimmomatic/trimmomatic_paper/analysis/
#$ -j y
#$ -o /mnt/data/project/2025_trimmomatic/trimmomatic_paper/logs/human_adapter_assessment.log
#$ -pe smp 40
#$ -l h_vmem=1G
#$ -l h_rt=02:00:00

# --- Environment Setup ---
set +e
set -u

source /mnt/data/personal/mambaforge/etc/profile.d/conda.sh
conda activate /mnt/data/project/2025_trimmomatic/trimmomatic_paper/scripts/parallel-gzip

# --- Configuration ---
ANALYSIS_DIR="/mnt/data/project/2025_trimmomatic/trimmomatic_paper/analysis/tools_human"
ADAPTER_SEED="AGATCGGAAGAGC"
THREADS=40

# --- Tool Detection ---
# Define Decompressor for .gz files
if command -v pragzip &> /dev/null; then
    DECOMPRESSor="pragzip -q -d -c -P $THREADS"
elif command -v pigz &> /dev/null; then
    DECOMPRESSor="pigz -dc -p $THREADS"
else
    DECOMPRESSor="zcat"
fi

# Define Searcher
if command -v rg &> /dev/null; then
    SEARCH="rg -c -F"
else
    SEARCH="grep -c -F"
fi

echo "--------------------------------------------------------------------------------"
echo "RESIDUAL ADAPTER ASSESSMENT"
echo "Searching for seed: $ADAPTER_SEED"
echo "--------------------------------------------------------------------------------"
printf "%-20s | %-15s | %-15s\n" "Tool" "R1 Matches" "R2 Matches"
echo "---------------------|-----------------|-----------------"

# --- File Definitions ---
declare -A FILES_R1
declare -A FILES_R2

# Trimmomatic
FILES_R1["Trimmomatic"]="${ANALYSIS_DIR}/trimmomatic_R1.fq.gz"
FILES_R2["Trimmomatic"]="${ANALYSIS_DIR}/trimmomatic_R2.fq.gz"

# fastp
FILES_R1["fastp"]="${ANALYSIS_DIR}/fastp_R1.fq.gz"
FILES_R2["fastp"]="${ANALYSIS_DIR}/fastp_R2.fq.gz"

# BBDuk
FILES_R1["BBDuk"]="${ANALYSIS_DIR}/bbduk_R1.fq.gz"
FILES_R2["BBDuk"]="${ANALYSIS_DIR}/bbduk_R2.fq.gz"

# Skewer
FILES_R1["Skewer"]="${ANALYSIS_DIR}/skewer_out-trimmed-pair1.fastq.gz"
FILES_R2["Skewer"]="${ANALYSIS_DIR}/skewer_out-trimmed-pair2.fastq.gz"

# Cutadapt
FILES_R1["Cutadapt"]="${ANALYSIS_DIR}/cutadapt_R1.fq.gz"
FILES_R2["Cutadapt"]="${ANALYSIS_DIR}/cutadapt_R2.fq.gz"

# RabbitTrim
FILES_R1["RabbitTrim"]="${ANALYSIS_DIR}/rabbit_out.read1_p"
FILES_R2["RabbitTrim"]="${ANALYSIS_DIR}/rabbit_out.read2_p"

# --- Execution Loop ---
for tool in "Trimmomatic" "fastp" "BBDuk" "Skewer" "Cutadapt" "RabbitTrim"; do
    unset count1 count2
    r1="${FILES_R1[$tool]}"
    r2="${FILES_R2[$tool]}"

    # Function to choose cat vs decompressor
    get_reader() {
        if [[ "$1" == *.gz ]]; then
            echo "$DECOMPRESSor"
        else
            echo "cat"
        fi
    }

    if [[ -f "$r1" ]]; then
        # Determine appropriate reader for R1 and R2
        READER1=$(get_reader "$r1")
        READER2=$(get_reader "$r2")

        count1=$($READER1 "$r1" | $SEARCH "$ADAPTER_SEED" 2>/dev/null || echo "0")
        count2=$($READER2 "$r2" | $SEARCH "$ADAPTER_SEED" 2>/dev/null || echo "0")
        
        printf "%-20s | %-15s | %-15s\n" "$tool" "$count1" "$count2"
    else
        printf "%-20s | %-15s | %-15s\n" "$tool" "MISSING FILE" "---"
    fi
done

echo "--------------------------------------------------------------------------------"