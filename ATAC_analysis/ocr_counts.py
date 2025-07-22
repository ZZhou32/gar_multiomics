import os
import glob

def extract_replicate_name(filename):
    """Extract replicate name (prefix before '_sorted_no_mito_')."""
    base = os.path.basename(filename)
    return base.split('_sorted_no_mito_')[0]

def count_lines_in_file(filepath):
    """Count non-empty lines in a file."""
    with open(filepath, 'r') as f:
        return sum(1 for line in f if line.strip())

def main():
    # Set your directory here (modify as needed)
    directory = "/mnt/research/FishEvoDevoGeno/raw_reads_4_tony/ATAC_peaks/ncOCR_narrow_peaks"  # Example: path to folder with .narrowPeak files
    
    # Find all .narrowPeak files in the specified directory (recursive)
    #narrowPeak_files = glob.glob(os.path.join(directory, '**/*.narrowPeak'), recursive=True)
    narrowPeak_files = glob.glob(os.path.join(directory, '**/*.bed'), recursive=True)

    if not narrowPeak_files:
        print(f"No .narrowPeak files found in: {directory}")
        return
    
    # Collect data: [(replicate_name, peak_count), ...]
    data = []
    for filepath in narrowPeak_files:
        replicate = extract_replicate_name(filepath)
        peak_count = count_lines_in_file(filepath)
        data.append((replicate, peak_count))
    
    # Write to TSV
    output_file = "/mnt/research/FishEvoDevoGeno/raw_reads_4_tony/ATAC_peaks/ncOCR_narrow_peaks/counts_ncOCR.tsv"
    with open(output_file, 'w') as f:
        f.write("Replicate\t# OCR\n")  # Header
        for replicate, count in data:
            f.write(f"{replicate}\t{count}\n")
    
    print(f"Results saved to: {output_file}")

if __name__ == "__main__":
    main()