#!/usr/bin/env python3
import sys
import os

# --- Configuration -----------------------------------------------------------
# The standard Illumina TruSeq 13bp Seed
SEED = "AGATCGGAAGAGC"

# EXPECTED EXTENSIONS (The bases immediately following the seed)
# R1 (Read 1): Reads into the Index Adapter (Starts with CACACGT...)
# We use a short 6bp signature to allow for some flexibility/errors
R1_EXTENSION_SIGS = ["CACACG", "CACGAC", "CACTCT", "CACC"] # Common variations

# R2 (Read 2): Reads into the Universal Adapter (Starts with GTCGTG...)
R2_EXTENSION_SIGS = ["GTCGTG", "GTCGGC", "GTCTCG"] 

# Minimum position: If adapter starts before base 15, assume it's a dimer/artifact
MIN_POS_THRESHOLD = 15 
# -----------------------------------------------------------------------------

def hamming_distance(s1, s2):
    """Calculate mismatches between two strings of equal length."""
    if len(s1) != len(s2): return 999
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))

def check_extension(sequence, seed_end_idx, is_read1):
    """
    Verifies if the bases following the seed match the expected adapter tail.
    Returns True if it looks like a real adapter.
    """
    # Extract the 6 bases immediately following the seed
    extension_obs = sequence[seed_end_idx : seed_end_idx + 6]
    
    # If read is too short to have an extension, it's ambiguous (treat as False or check context)
    if len(extension_obs) < 4:
        return False

    valid_sigs = R1_EXTENSION_SIGS if is_read1 else R2_EXTENSION_SIGS
    
    # Check against valid signatures (Allow 1 mismatch for sequencing error)
    for sig in valid_sigs:
        # Truncate sig if extension is shorter (e.g. adapter at very end of read)
        comp_len = min(len(sig), len(extension_obs))
        if hamming_distance(extension_obs[:comp_len], sig[:comp_len]) <= 1:
            return True
            
    return False

def analyze_file(filepath, is_read1):
    total_hits = 0
    true_residuals = 0
    fp_extension = 0
    fp_position = 0
    
    tool_name = os.path.basename(filepath).split('.')[0]
    read_type = "R1" if is_read1 else "R2"

    with open(filepath, 'r') as f:
        lines = f.readlines()

    # Process in chunks of 3 lines (Header, Sequence, Separator --)
    # Based on your grep output format
    for i in range(0, len(lines)):
        line = lines[i].strip()
        
        # We only care about the sequence line (not header starting with @, not separator --)
        if line.startswith("@") or line.startswith("--") or len(line) < 20:
            continue
            
        # Find the seed in this sequence
        # Note: grep might output multiple matches, we take the first valid one
        idx = line.find(SEED)
        
        if idx != -1:
            total_hits += 1
            
            # Check 1: Position (Is it a dimer/short fragment?)
            if idx < MIN_POS_THRESHOLD:
                fp_position += 1
                continue
            
            # Check 2: Extension (Is it the correct adapter?)
            seed_end = idx + len(SEED)
            if check_extension(line, seed_end, is_read1):
                true_residuals += 1
                # Optional: Print the true positive for verification
                # print(f"[TP] {line[idx:]}") 
            else:
                fp_extension += 1
                # Optional: Print the mismatch for debugging
                # print(f"[FP] {line[idx:idx+20]}... Expected {R1_EXTENSION_SIGS if is_read1 else R2_EXTENSION_SIGS}")

    return {
        "Tool": tool_name,
        "Read": read_type,
        "Total_Candidates": total_hits,
        "True_Residuals": true_residuals,
        "FP_Random_Match": fp_extension,
        "FP_Dimer_Start": fp_position
    }

# --- Main Execution ----------------------------------------------------------
if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python3 verify_residuals.py <R1_file> <R2_file> [More files...]")
        sys.exit(1)

    print(f"{'Tool':<15} | {'Read':<4} | {'Candidates':<10} | {'TRUE RESIDUALS':<15} | {'False Positives':<15}")
    print("-" * 75)

    for fpath in sys.argv[1:]:
        if not os.path.exists(fpath): continue
        
        # Auto-detect R1 vs R2 from filename
        is_r1 = "R1" in fpath or "_1" in fpath
        
        stats = analyze_file(fpath, is_r1)
        
        print(f"{stats['Tool']:<15} | {stats['Read']:<4} | "
              f"{stats['Total_Candidates']:<10} | "
              f"{stats['True_Residuals']:<15} | "
              f"{stats['FP_Random_Match'] + stats['FP_Dimer_Start']:<15}")