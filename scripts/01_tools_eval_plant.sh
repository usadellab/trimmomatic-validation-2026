#!/bin/bash

#-------------------------------------------------------------------------------
# SGE Job Script: Trimming Tools Benchmarking (2026 Revision)
#
# TOOLS EVALUATED:
# 1. Trimmomatic v0.40 (Java)
# 2. fastp (C++ / Conda)
# 3. BBDuk (Java / Conda)
# 4. Skewer (C++ / Conda)
# 5. Cutadapt (Python / Conda)
# 6. RabbitTrim (C++ / Native)
#-------------------------------------------------------------------------------

#$ -N tools_eval
#$ -wd /mnt/data/project/2025_trimmomatic/trimmomatic_paper/analysis/
#$ -j y
#$ -o /mnt/data/project/2025_trimmomatic/trimmomatic_paper/logs/tools_eval.log
#$ -pe smp 40
#$ -l h_vmem=3.125G
#$ -l h_rt=24:00:00

# --- Environment Setup ---
set +e              # Disable exit-on-error
set +o pipefail     # Disable pipefail
set -u              # Treat unset variables as an error

# Raise system limits
ulimit -u 65535 2>/dev/null || true
ulimit -n 65535 2>/dev/null || true

source /mnt/data/personal/mambaforge/etc/profile.d/conda.sh

# --- Configuration ---
LOG_DIR="/mnt/data/project/2025_trimmomatic/trimmomatic_paper/logs"
RAW_DIR="/mnt/data/project/2025_trimmomatic/trimmomatic_paper/raw"
ANALYSIS_DIR="/mnt/data/project/2025_trimmomatic/trimmomatic_paper/analysis/tools"
mkdir -p "$LOG_DIR" "$ANALYSIS_DIR"

# Input Files
READ1="${RAW_DIR}/MS01_CasR09_R1_001.fastq.gz"
READ2="${RAW_DIR}/MS01_CasR09_R2_001.fastq.gz"
ADAPTERS_FILE="/mnt/bin/trimmomatic/Trimmomatic-0.40/adapters/TruSeq3-PE-2-GGGGG.fa"
TRIMMOMATIC_JAR="/mnt/bin/trimmomatic/Trimmomatic-0.40/trimmomatic-0.40.jar"

# Parameters
THREADS=40
TRIM_STEPS="ILLUMINACLIP:${ADAPTERS_FILE}:2:30:10 SLIDINGWINDOW:4:20 MINLEN:36"

# Check for skip argument
SKIP_EXISTING=false
if [[ "${1:-}" == "--skip-existing" ]]; then
    SKIP_EXISTING=true
fi

log() {
  echo "[$(date +'%Y-%m-%d %H:%M:%S')] $1"
}

# Logic: Skip only if flag is set AND log file has content (>0 bytes)
should_skip() {
    local log_file="$1"
    if [[ "$SKIP_EXISTING" == "true" ]] && [[ -s "$log_file" ]]; then
        return 0 # True, skip
    else
        return 1 # False, run
    fi
}

log "Starting Benchmarking on Host: $(hostname)"
log "System Limits: Max Processes=$(ulimit -u)"

# --- 1. Trimmomatic v0.40 ---
LOG_FILE="${LOG_DIR}/trimmomatic.time.log"
if should_skip "$LOG_FILE"; then
    log "Skipping Trimmomatic (Log exists)."
else
    log "Running Trimmomatic v0.40..."
    # Defined outputs for Paired and Unpaired streams
    OUT_R1_P="${ANALYSIS_DIR}/trimmomatic_R1.fq.gz"
    OUT_R1_U="${ANALYSIS_DIR}/trimmomatic_R1_unpaired.fq.gz"
    OUT_R2_P="${ANALYSIS_DIR}/trimmomatic_R2.fq.gz"
    OUT_R2_U="${ANALYSIS_DIR}/trimmomatic_R2_unpaired.fq.gz"

    /usr/bin/time -v java -jar "$TRIMMOMATIC_JAR" PE -threads "$THREADS" \
        "$READ1" "$READ2" \
        "$OUT_R1_P" "$OUT_R1_U" \
        "$OUT_R2_P" "$OUT_R2_U" \
        $TRIM_STEPS \
        &> "$LOG_FILE"

    if [[ $? -eq 0 ]]; then log "Trimmomatic success."; else log "ERROR: Trimmomatic failed."; fi
fi

# --- 2. fastp ---
LOG_FILE="${LOG_DIR}/fastp.time.log"
if should_skip "$LOG_FILE"; then
    log "Skipping fastp (Log exists)."
else
    log "Running fastp..."
    if conda activate /mnt/bin/fastp/fastp-v1.1.0-conda; then
        /usr/bin/time -v fastp --thread "$THREADS" \
            -i "$READ1" -I "$READ2" \
            -o "${ANALYSIS_DIR}/fastp_R1.fq.gz" -O "${ANALYSIS_DIR}/fastp_R2.fq.gz" \
            --adapter_fasta "$ADAPTERS_FILE" \
            -j "${LOG_DIR}/fastp.json" -h "${LOG_DIR}/fastp.html" \
            &> "$LOG_FILE"
        
        if [[ $? -eq 0 ]]; then log "fastp success."; else log "ERROR: fastp failed."; fi
        conda deactivate
    else
        log "ERROR: Could not activate fastp environment."
    fi
fi

# --- 3. BBDuk ---
LOG_FILE="${LOG_DIR}/bbduk.time.log"
if should_skip "$LOG_FILE"; then
    log "Skipping BBDuk (Log exists)."
else
    log "Running BBDuk..."
    if conda activate /mnt/bin/bbmap/bbmap-v39.52-conda; then
        /usr/bin/time -v bbduk.sh -Xmx90g threads="$THREADS" ordered=f \
            in="$READ1" in2="$READ2" \
            out="${ANALYSIS_DIR}/bbduk_R1.fq.gz" out2="${ANALYSIS_DIR}/bbduk_R2.fq.gz" \
            ref="$ADAPTERS_FILE" ktrim=r k=23 mink=11 hdist=1 tpe tbo \
            &> "$LOG_FILE"

        if [[ $? -eq 0 ]]; then log "BBDuk success."; else log "ERROR: BBDuk failed."; fi
        conda deactivate
    else
        log "ERROR: Could not activate BBDuk Conda environment."
    fi
fi

# --- 4. Skewer ---
LOG_FILE="${LOG_DIR}/skewer.time.log"
if should_skip "$LOG_FILE"; then
    log "Skipping Skewer (Log exists)."
else
    log "Running Skewer..."
    if conda activate /mnt/bin/skewer/skewer-v0.2.2-conda; then
        /usr/bin/time -v skewer -t "$THREADS" -m pe -z -x "$ADAPTERS_FILE" \
            -o "${ANALYSIS_DIR}/skewer_out" "$READ1" "$READ2" \
            &> "$LOG_FILE"

        if [[ $? -eq 0 ]]; then log "Skewer success."; else log "ERROR: Skewer failed."; fi
        conda deactivate
    else
        log "ERROR: Could not activate Skewer environment."
    fi
fi

# --- 5. Cutadapt ---
LOG_FILE="${LOG_DIR}/cutadapt.time.log"
if should_skip "$LOG_FILE"; then
    log "Skipping Cutadapt (Log exists)."
else
    log "Running Cutadapt..."
    if conda activate /mnt/bin/cutadapt/cutadapt-v5.2-conda; then
        /usr/bin/time -v cutadapt -j "$THREADS" \
            -a "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC" \
            -a "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG" \
            -A "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA" \
            -A "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG" \
            -o "${ANALYSIS_DIR}/cutadapt_R1.fq.gz" -p "${ANALYSIS_DIR}/cutadapt_R2.fq.gz" \
            "$READ1" "$READ2" \
            &> "$LOG_FILE"

        if [[ $? -eq 0 ]]; then log "Cutadapt success."; else log "ERROR: Cutadapt failed."; fi
        conda deactivate
    else
        log "ERROR: Could not activate Cutadapt environment."
    fi
fi

# --- 6. RabbitTrim ---
LOG_FILE="${LOG_DIR}/rabbittrim.time.log"
if should_skip "$LOG_FILE"; then
    log "Skipping RabbitTrim (Log exists)."
else
    log "Running RabbitTrim..."    
    /usr/bin/time -v /mnt/bin/rabbittrim/RabbitTrim_v2.0.0/build/RabbitTrim trimmomatic \
        --PE \
        --threads "$THREADS" \
        --forward "$READ1" \
        --reverse "$READ2" \
        --output "${ANALYSIS_DIR}/rabbit_out" \
        --stats "${LOG_DIR}/rabbittrim.stats" \
        --steps $TRIM_STEPS \
        --compressLevel 1 \
        &> "$LOG_FILE"

    if [[ $? -eq 0 ]]; then log "RabbitTrim success."; else log "ERROR: RabbitTrim failed."; fi
fi

log "All benchmarks complete."