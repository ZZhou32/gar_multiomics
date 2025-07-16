#!/usr/bin/env python3
import os
import subprocess
from pathlib import Path

# ===== Configuration =====
# Input directory containing MACS2 subfolders
MACS2_DIR = "/mnt/research/FishEvoDevoGeno/raw_reads_4_tony/ATAC_bam/MACS2"

# Output directories
RAW_PEAKS_DIR = "/mnt/research/FishEvoDevoGeno/raw_reads_4_tony/ATAC_peaks/narrowPeaks_raw"
CLEAN_PEAKS_DIR = "/mnt/research/FishEvoDevoGeno/raw_reads_4_tony/ATAC_peaks/ncOCR_narrow_peaks"

# Exclusion regions (must be BED format)
EXCLUSION_BED = "/mnt/research/FishEvoDevoGeno/raw_reads_4_tony/exons_only.bed"

# Bedtools subtract options
BEDTOOLS_OPTS = "-A"  # Remove entire peak if any overlap

# ===== Main Script =====
def main():
    # Create output directories if they don't exist
    os.makedirs(RAW_PEAKS_DIR, exist_ok=True)
    os.makedirs(CLEAN_PEAKS_DIR, exist_ok=True)
    
    print(f"1. Collecting narrowPeak files from {MACS2_DIR}...")
    narrowpeak_files = list(Path(MACS2_DIR).rglob("*.narrowPeak"))
    
    if not narrowpeak_files:
        raise FileNotFoundError(f"No .narrowPeak files found in {MACS2_DIR}")
    
    print(f"2. Copying {len(narrowpeak_files)} files to {RAW_PEAKS_DIR}...")
    for src_file in narrowpeak_files:
        dest_file = Path(RAW_PEAKS_DIR) / src_file.name
        if not dest_file.exists():
            subprocess.run(["cp", str(src_file), str(dest_file)], check=True)
    
    print("3. Processing files with bedtools subtract...")
    for raw_file in Path(RAW_PEAKS_DIR).glob("*.narrowPeak"):
        output_file = Path(CLEAN_PEAKS_DIR) / f"{raw_file.stem}_ncOCR.bed"
        
        cmd = f"bedtools subtract {BEDTOOLS_OPTS} -a {raw_file} -b {EXCLUSION_BED} > {output_file}"
        
        print(f"Processing {raw_file.name} → {output_file.name}")
        try:
            subprocess.run(cmd, shell=True, check=True, executable="/bin/bash")
        except subprocess.CalledProcessError as e:
            print(f"Error processing {raw_file}: {e}")
            continue
    
    print(f"\nDone! Cleaned peaks saved to {CLEAN_PEAKS_DIR}")

if __name__ == "__main__":
    main()