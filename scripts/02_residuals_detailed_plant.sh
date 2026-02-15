#!/bin/bash
#$ -N adapter_assessment_detailed
#$ -wd /mnt/data/project/2025_trimmomatic/trimmomatic_paper/analysis/
#$ -j y
#$ -o /mnt/data/project/2025_trimmomatic/trimmomatic_paper/logs/adapter_assessment_detailed_plant.log
#$ -pe smp 40
#$ -l h_vmem=1G
#$ -l h_rt=02:00:00

# --- Environment Setup ---
set +e
set -u

source /mnt/data/personal/mambaforge/etc/profile.d/conda.sh
# Ensure your environment has pigz/pragzip and rg (ripgrep) if possible
conda activate /mnt/data/project/2025_trimmomatic/trimmomatic_paper/scripts/parallel-gzip

# --- Configuration ---
ANALYSIS_DIR="/mnt/data/project/2025_trimmomatic/trimmomatic_paper/analysis/scaling"
RESIDUAL_OUT_DIR="/mnt/data/project/2025_trimmomatic/trimmomatic_paper/analysis/residuals_detailed_plant"
ADAPTER_SEED="AGATCGGAAGAGC"
THREADS=40

mkdir -p "$RESIDUAL_OUT_DIR"

# --- Tool Detection ---
# 1. Decompressor
if command -v pragzip &> /dev/null; then
    DECOMPRESS="pragzip -q -d -c -P $THREADS"
elif command -v pigz &> /dev/null; then
    DECOMPRESS="pigz -dc -p $THREADS"
else
    DECOMPRESS="zcat"
fi

# 2. Searcher (Must support -B 1 to get Read ID)
# We prefer ripgrep (rg) for speed, falling back to grep
if command -v rg &> /dev/null; then
    # -B 1: Context Before (Header), -F: Fixed string, -N: No line numbers
    SEARCH_CMD="rg -B 1 -F --no-line-number --no-heading"
else
    SEARCH_CMD="grep -B 1 -F"
fi

echo "--------------------------------------------------------------------------------"
echo "RESIDUAL ADAPTER ASSESSMENT (DETAILED)"
echo "Searching for seed: $ADAPTER_SEED"
echo "Output Directory: $RESIDUAL_OUT_DIR"
echo "--------------------------------------------------------------------------------"
printf "%-20s | %-15s | %-15s\n" "Tool" "R1 Matches" "R2 Matches"
echo "---------------------|-----------------|-----------------"

# --- File Definitions ---
declare -A FILES_R1
declare -A FILES_R2

# Trimmomatic
FILES_R1["Trimmomatic"]="${ANALYSIS_DIR}/Plant_Trimmomatic_40t_R1.fq.gz"
FILES_R2["Trimmomatic"]="${ANALYSIS_DIR}/Plant_Trimmomatic_40t_R2.fq.gz"

# fastp
FILES_R1["fastp"]="${ANALYSIS_DIR}/Plant_fastp_40t_R1.fq.gz"
FILES_R2["fastp"]="${ANALYSIS_DIR}/Plant_fastp_40t_R2.fq.gz"

# BBDuk
FILES_R1["BBDuk"]="${ANALYSIS_DIR}/Plant_BBDuk_40t_R1.fq.gz"
FILES_R2["BBDuk"]="${ANALYSIS_DIR}/Plant_BBDuk_40t_R2.fq.gz"

# Skewer
FILES_R1["Skewer"]="${ANALYSIS_DIR}/Plant_skewer_40t_out-trimmed-pair1.fastq.gz"
FILES_R2["Skewer"]="${ANALYSIS_DIR}/Plant_skewer_40t_out-trimmed-pair2.fastq.gz"

# Cutadapt
FILES_R1["Cutadapt"]="${ANALYSIS_DIR}/Plant_Cutadapt_40t_R1.fq.gz"
FILES_R2["Cutadapt"]="${ANALYSIS_DIR}/Plant_Cutadapt_40t_R2.fq.gz"

# RabbitTrim
FILES_R1["RabbitTrim"]="${ANALYSIS_DIR}/Plant_RabbitTrim_40t_out.read1_p.gz"
FILES_R2["RabbitTrim"]="${ANALYSIS_DIR}/Plant_RabbitTrim_40t_out.read2_p.gz"

# --- Execution Loop ---
for tool in "Trimmomatic" "fastp" "BBDuk" "Skewer" "Cutadapt" "RabbitTrim"; do
    unset count1 count2
    r1="${FILES_R1[$tool]}"
    r2="${FILES_R2[$tool]}"
    
    # Define Output Files for Residuals
    out_r1="${RESIDUAL_OUT_DIR}/${tool}.R1.residuals.txt"
    out_r2="${RESIDUAL_OUT_DIR}/${tool}.R2.residuals.txt"

    # Helper for cat vs unzip
    get_reader() {
        if [[ "$1" == *.gz ]]; then echo "$DECOMPRESS"; else echo "cat"; fi
    }

    if [[ -f "$r1" ]]; then
        READER1=$(get_reader "$r1")
        READER2=$(get_reader "$r2")

        # 1. Search R1
        # Run search, save to file
        $READER1 "$r1" | $SEARCH_CMD "$ADAPTER_SEED" > "$out_r1"
        # Count occurrences (grep -c is faster than wc -l / 2)
        # We count the seed itself in the output file to avoid counting headers
        count1=$(grep -c "$ADAPTER_SEED" "$out_r1")

        # 2. Search R2
        $READER2 "$r2" | $SEARCH_CMD "$ADAPTER_SEED" > "$out_r2"
        count2=$(grep -c "$ADAPTER_SEED" "$out_r2")
        
        printf "%-20s | %-15s | %-15s\n" "$tool" "$count1" "$count2"
    else
        printf "%-20s | %-15s | %-15s\n" "$tool" "MISSING FILE" "---"
    fi
done

echo "--------------------------------------------------------------------------------"
echo "Detailed residual files saved to: $RESIDUAL_OUT_DIR"