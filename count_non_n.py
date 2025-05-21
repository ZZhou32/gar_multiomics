import gzip
import sys

def count_non_n_bases(fasta_path):
    """
    Counts all non-N bases in a FASTA file (gzipped or plain text)
    """
    total = 0
    open_func = gzip.open if fasta_path.endswith('.gz') else open
    
    try:
        with open_func(fasta_path, 'rt') as f:  # 'rt' mode for text reading
            for line in f:
                if not line.startswith('>'):  # Skip header lines
                    total += len(line.strip()) - line.upper().count('N')
    except Exception as e:
        print(f"Error reading file: {e}", file=sys.stderr)
        return None
        
    return total

if __name__ == "__main__":
    fasta_file = "/mnt/research/FishEvoDevoGeno/raw_reads_4_tony/GCF_040954835.1_fLepOcu1.hap2_genomic.fna"
    
    print(f"Counting non-N bases in {fasta_file}...")
    count = count_non_n_bases(fasta_file)
    
    if count is not None:
        print(f"Total non-N bases: {count:,}")
        print(f"Approximate genome size (bp): {count:,}")
    else:
        print("Failed to count bases", file=sys.stderr)