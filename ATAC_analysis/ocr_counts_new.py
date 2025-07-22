import os
import glob
import subprocess
import re
from collections import defaultdict
### change def main based on ur need.

def extract_stage_and_replicate(filename):
    """Extract stage and replicate information from filename."""
    base = os.path.basename(filename)
    
    # Handle the complex filename pattern
    # Example: "Gar_ATAC_whole_embryo_32-33_rep3_sorted_no_mito_resorted_no_dup_uniq_query..."
    # Should extract:
    # - Stage: Gar_ATAC_whole_embryo_32-33
    # - Replicate: rep3
    
    # First remove everything after first occurrence of '_sorted_no_mito_'
    if '_sorted_no_mito_' in base:
        base = base.split('_sorted_no_mito_')[0]
    
    # Now look for _rep\d+ pattern
    rep_match = re.search(r'_rep(\d+)', base)
    if rep_match:
        # Split at the _rep\d+ pattern
        stage = base[:rep_match.start()]
        replicate = f"rep{rep_match.group(1)}"
    else:
        # No replicate number found
        stage = base
        replicate = "rep1"  # Default for single replicates
    
    return stage, replicate

def count_lines_in_file(filepath):
    """Count non-empty lines in a file."""
    with open(filepath, 'r') as f:
        return sum(1 for line in f if line.strip())

def merge_peaks_with_bedtools(input_files, output_file):
    """Merge peaks using bedtools merge with proper sorting."""
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    
    # 1. Sort each file with bedtools (handles chr names correctly)
    sorted_files = []
    for f in input_files:
        sorted_file = f"{f}.sorted"
        cmd_sort = f"bedtools sort -i {f} > {sorted_file}"
        subprocess.run(cmd_sort, shell=True, check=True)
        sorted_files.append(sorted_file)
    
    # 2. Concatenate all sorted files
    cat_file = f"{output_file}.cat"
    cmd_cat = f"cat {' '.join(sorted_files)} > {cat_file}"
    subprocess.run(cmd_cat, shell=True, check=True)
    
    # 3. Sort the concatenated file (final safety check)
    cat_sorted = f"{cat_file}.sorted"
    cmd_sort_cat = f"bedtools sort -i {cat_file} > {cat_sorted}"
    subprocess.run(cmd_sort_cat, shell=True, check=True)
    
    # 4. Merge with bedtools
    cmd_merge = f"bedtools merge -i {cat_sorted} > {output_file}"
    subprocess.run(cmd_merge, shell=True, check=True)
    
    # Cleanup
    for f in sorted_files + [cat_file, cat_sorted]:
        try:
            os.remove(f)
        except FileNotFoundError:
            pass

def generate_stage_specific_counts(output_dir):
    """Generate OCR counts for stage-specific merged peaks."""
    # Find all merged peak files in the stage-specific directory
    merged_peak_files = glob.glob(os.path.join(output_dir, '*_combined_peaks.bed'))
    
    if not merged_peak_files:
        print(f"No merged peak files found in: {output_dir}")
        return
    
    # Create counts file in the stage-specific directory
    counts_file = os.path.join(output_dir, "stage_specific_OCR_counts.txt")
    
    with open(counts_file, 'w') as f:
        f.write("Stage\t# OCR (merged peaks)\tFile\n")  # Header
        
        for filepath in merged_peak_files:
            stage = os.path.basename(filepath).replace('_combined_peaks.bed', '')
            peak_count = count_lines_in_file(filepath)
            f.write(f"{stage}\t{peak_count}\t{os.path.basename(filepath)}\n")
    
    print(f"\nStage-specific OCR counts saved to: {counts_file}")

def main():
    # Set your directory here (modify as needed)
    directory = "/mnt/research/FishEvoDevoGeno/raw_reads_4_tony/ATAC_peaks/ncOCR_narrow_peaks"
    output_dir = os.path.join(directory, "stage-specific")
    
    # Find all .bed files in the specified directory (recursive)
    peak_files = glob.glob(os.path.join(directory, '**/*.bed'), recursive=True)
    #peak_files = glob.glob(os.path.join(directory, '**/*.narrowPeak'), recursive=True)
    if not peak_files:
        print(f"No .bed files found in: {directory}")
        return
    
    # Group files by stage and collect data
    stage_files = defaultdict(list)
    data = []
    
    for filepath in peak_files:
        stage, replicate = extract_stage_and_replicate(filepath)
        stage_files[stage].append(filepath)
        peak_count = count_lines_in_file(filepath)
        data.append((stage, replicate, peak_count, filepath))
    
    # Print some debug info
    print("Found the following stage/replicate combinations:")
    for stage, replicate, count, path in data:
        print(f"- Stage: {stage}, Replicate: {replicate}, Count: {count}, File: {os.path.basename(path)}")
    
    # Write counts to TSV
    counts_file = os.path.join(directory, "counts_ncOCR.tsv")
    with open(counts_file, 'w') as f:
        f.write("Stage\tReplicate\t# OCR\tFile\n")  # Header
        for stage, replicate, count, path in data:
            f.write(f"{stage}\t{replicate}\t{count}\t{os.path.basename(path)}\n")
    
    print(f"\nCount results saved to: {counts_file}")
    
    # Merge peaks for each stage
    print("\nMerging peaks by stage:")
    for stage, files in stage_files.items():
        output_file = os.path.join(output_dir, f"{stage}_combined_peaks.bed")
        print(f"\nProcessing {len(files)} files for stage '{stage}':")
        for f in files:
            print(f"- {os.path.basename(f)}")
        
        merge_peaks_with_bedtools(files, output_file)
        print(f"Saved merged peaks to: {output_file}")
    
    # Generate stage-specific OCR counts
    generate_stage_specific_counts(output_dir)
    
    print("\nAll stages processed.")

if __name__ == "__main__":
    main()