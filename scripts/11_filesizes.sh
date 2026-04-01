#!/bin/bash

# --- Configuration ---
OUT_DIR="/mnt/data/project/2025_trimmomatic/trimmomatic_paper/analysis/scaling_extended_v2"
DATASETS=("Plant" "Human")
TOOLS=("Trimmomatic" "RabbitTrim" "fastp" "BBDuk" "skewer" "Cutadapt")
THREAD_COUNTS=(8 40)

echo -e "Tool\tDataset\tThreads\tRead\tSize_Bytes\tSize_MB"

for d in "${DATASETS[@]}"; do
    for t in "${TOOLS[@]}"; do
        for tc in "${THREAD_COUNTS[@]}"; do
            
            case "$t" in
                "RabbitTrim")
                    r1_file="${OUT_DIR}/${d}_${t}_${tc}t_out.read1_p.gz"
                    r2_file="${OUT_DIR}/${d}_${t}_${tc}t_out.read2_p.gz"
                    ;;
                "skewer")
                    r1_file="${OUT_DIR}/${d}_${t}_${tc}t_out-trimmed-pair1.fastq.gz"
                    r2_file="${OUT_DIR}/${d}_${t}_${tc}t_out-trimmed-pair2.fastq.gz"
                    ;;
                *)
                    r1_file="${OUT_DIR}/${d}_${t}_${tc}t_R1.fq.gz"
                    r2_file="${OUT_DIR}/${d}_${t}_${tc}t_R2.fq.gz"
                    ;;
            esac

            # Log size for Read 1
            if [[ -f "$r1_file" ]]; then
                size1=$(stat -c%s "$r1_file")
                size1_mb=$(echo "scale=2; $size1/1048576" | bc)
                echo -e "${t}\t${d}\t${tc}\tR1\t${size1}\t${size1_mb}"
            fi

            # Log size for Read 2
            if [[ -f "$r2_file" ]]; then
                size2=$(stat -c%s "$r2_file")
                size2_mb=$(echo "scale=2; $size2/1048576" | bc)
                echo -e "${t}\t${d}\t${tc}\tR2\t${size2}\t${size2_mb}"
            fi
        done
    done
done