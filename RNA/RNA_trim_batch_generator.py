import os

def find_paired_files_and_bases(folder_path):
    """
    Scan folder for all .fq.gz files, find their pairs (_1 and _2), 
    and extract base names for each complete pair.
    """
    if not os.path.exists(folder_path):
        print(f"Error: Folder {folder_path} does not exist!")
        return {}
    
    # Get all .fq.gz files
    all_files = [f for f in os.listdir(folder_path) if f.endswith('.fq.gz')]
    print(f"Found {len(all_files)} .fq.gz files in total")
    
    # Dictionary to store file pairs: {base_name: {'R1': file1, 'R2': file2}}
    paired_files = {}
    
    # Process each file to find pairs
    for file_name in all_files:
        print(f"Processing: {file_name}")
        
        # Check if this is an R1 file (ends with _1.fq.gz)
        if file_name.endswith('_1.fq.gz'):
            # Extract base name by removing _1.fq.gz
            base_name = file_name[:-7]  # Remove '_1.fq.gz'
            
            # Look for corresponding R2 file
            r2_file = base_name + '2.fq.gz'
            
            if r2_file in all_files:
                paired_files[base_name] = {
                    'R1': file_name,
                    'R2': r2_file,
                    'R1_path': os.path.join(folder_path, file_name),
                    'R2_path': os.path.join(folder_path, r2_file)
                }
                print(f"  ✓ Found pair: {base_name}")
                print(f"    R1: {file_name}")
                print(f"    R2: {r2_file}")
            else:
                print(f"  ✗ No R2 pair found for {file_name} (looking for {r2_file})")
    
    print(f"\nSummary:")
    print(f"Total files scanned: {len(all_files)}")
    print(f"Complete pairs found: {len(paired_files)}")
    print(f"Unpaired files: {len(all_files) - (len(paired_files) * 2)}")
    
    if paired_files:
        print(f"\nPaired samples found:")
        for base_name in sorted(paired_files.keys()):
            print(f"  - {base_name}")
    
    return paired_files

def generate_rna_batch_files_and_master(folder_path, batch_folder, trim_folder):
    """
    Generate individual batch files for RNA-seq trimming based on paired files found.
    """
    # Find all paired files
    paired_files = find_paired_files_and_bases(folder_path)
    
    if not paired_files:
        print("No paired files found! Cannot generate batch scripts.")
        return

    # Ensure the output folders exist
    os.makedirs(batch_folder, exist_ok=True)
    os.makedirs(trim_folder, exist_ok=True)
    os.makedirs(f"{trim_folder}/out", exist_ok=True)
    print(f"\nCreated output directories:")
    print(f"  - Batch files: {batch_folder}")
    print(f"  - Trimmed output: {trim_folder}")
    print(f"  - Logs: {trim_folder}/out")

    # Open the master batch file for writing
    master_batch_path = os.path.join(batch_folder, "master_rna_trim.sb")
    
    with open(master_batch_path, "w") as master_batch_file:
        master_batch_file.write("#!/bin/bash --login\n\n")
        master_batch_file.write("# Master script to submit all RNA-seq trimming jobs\n")
        master_batch_file.write(f"# Generated for {len(paired_files)} paired samples\n")
        master_batch_file.write(f"# Source folder: {folder_path}\n\n")

        batch_files_created = 0
        
        for base_name in sorted(paired_files.keys()):
            file_info = paired_files[base_name]
            
            # Create batch file name using base name
            batch_file_name = os.path.join(batch_folder, f"{base_name}_trim.sb")
            
            print(f"Creating batch file for: {base_name}")

            # Customize the content of the batch file to match working script
            batch_content = f"""#!/bin/bash --login

#SBATCH --job-name={base_name}_fastp_job
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --mail-user=zhouzeh2@msu.edu
#SBATCH --output={trim_folder}/out/%x-%j.SLURMout
#SBATCH --error={trim_folder}/out/%x-%j.SLURMerr

echo "Processing: {base_name}"
echo "Job started at: $(date)"

# File paths
in1="{file_info['R1_path']}"
in2="{file_info['R2_path']}"
out1="{trim_folder}/{base_name}_1NOADAPTERR.fq.gz"
out2="{trim_folder}/{base_name}_2NOADAPTERR.fq.gz"
html_report="{trim_folder}/{base_name}_fastp_report.html"
json_report="{trim_folder}/{base_name}_fastp_report.json"

# Create output directory
mkdir -p {trim_folder}

echo "Input file sizes:"
ls -lh "$in1" "$in2"

# Set up fastp environment (custom installation)
export PATH=/mnt/home/zhouzeh2/software/bin:$PATH
export LD_LIBRARY_PATH=/mnt/home/zhouzeh2/software/lib:$LD_LIBRARY_PATH

echo "Memory before processing:"
free -h

echo "Using fastp version:"
fastp --version

# Verify input files exist
echo "Verifying input files..."
if [[ ! -f "$in1" ]]; then
    echo "ERROR: Input R1 file not found: $in1"
    exit 1
fi

if [[ ! -f "$in2" ]]; then
    echo "ERROR: Input R2 file not found: $in2"
    exit 1
fi

echo "✓ Both input files verified"

# Run fastp with minimal, safe settings (matching working script)
fastp \\
    --in1 "$in1" \\
    --in2 "$in2" \\
    --out1 "$out1" \\
    --out2 "$out2" \\
    --html "$html_report" \\
    --json "$json_report" \\
    --detect_adapter_for_pe \\
    --cut_tail \\
    --cut_window_size=4 \\
    --cut_mean_quality=20 \\
    --qualified_quality_phred=20 \\
    --length_required=36 \\
    --thread=1 \\
    --compression=1

if [[ $? -eq 0 ]]; then
    echo "✓ fastp completed successfully"
    echo "Output files:"
    ls -lh "$out1" "$out2"
    echo "Reports:"
    ls -lh "$html_report" "$json_report"
    
    # Get read counts
    INPUT_READS=$(zcat "$in1" | wc -l | awk '{{print int($1/4)}}')
    OUTPUT_READS=$(zcat "$out1" | wc -l | awk '{{print int($1/4)}}')
    echo "Input reads: $INPUT_READS"
    echo "Output reads: $OUTPUT_READS"
    
    if [[ $INPUT_READS -gt 0 ]]; then
        RETENTION=$(echo "scale=1; $OUTPUT_READS * 100 / $INPUT_READS" | bc)
        echo "Retention rate: $RETENTION%"
    fi
else
    echo "✗ fastp failed"
    exit 1
fi

echo "Job completed: $(date)"
"""
            
            # Write the batch file
            with open(batch_file_name, "w") as batch_file:
                batch_file.write(batch_content)
            
            # Make the batch file executable
            os.chmod(batch_file_name, 0o755)
            batch_files_created += 1

            # Add sbatch command to the master batch file
            master_batch_file.write(f"echo 'Submitting {base_name} RNA-seq trimming job...'\n")
            master_batch_file.write(f"sbatch {batch_file_name}\n")
            master_batch_file.write(f"echo 'Job submitted for {base_name}'\n")
            master_batch_file.write(f"sleep 2  # Delay between submissions\n\n")

        # Add summary to master batch file
        master_batch_file.write(f"echo 'All {batch_files_created} RNA-seq trimming jobs submitted!'\n")
        master_batch_file.write(f"echo 'Monitor with: squeue -u $USER'\n")
        master_batch_file.write(f"echo 'Check logs in: {trim_folder}/out/'\n")

    # Make master batch file executable
    os.chmod(master_batch_path, 0o755)
    
    print(f"\n" + "="*60)
    print(f"BATCH FILE GENERATION COMPLETED")
    print(f"="*60)
    print(f"Paired samples processed: {len(paired_files)}")
    print(f"Batch files created: {batch_files_created}")
    print(f"Master script: {master_batch_path}")
    print(f"Individual batch files: {batch_folder}")
    print(f"\nTo submit all jobs:")
    print(f"bash {master_batch_path}")
    print(f"\nTo submit individual jobs:")
    print(f"sbatch {batch_folder}/<sample_name>_trim.sb")

def preview_batch_generation(folder_path):
    """
    Preview what batch files would be generated without actually creating them.
    """
    print("="*60)
    print("PREVIEW MODE - No files will be created")
    print("="*60)
    
    paired_files = find_paired_files_and_bases(folder_path)
    
    if paired_files:
        print(f"\nWould generate {len(paired_files)} batch files:")
        for base_name in sorted(paired_files.keys()):
            file_info = paired_files[base_name]
            print(f"  {base_name}_trim.sb")
            print(f"    R1: {file_info['R1']}")
            print(f"    R2: {file_info['R2']}")
            print()
    
    return len(paired_files)

# Configuration
folder_path = "/mnt/research/FishEvoDevoGeno/raw_reads_4_tony/RNA" 
batch_output_folder = "/mnt/research/FishEvoDevoGeno/raw_reads_4_tony/script/RNA/trim" 
trim_output_folder = "/mnt/research/FishEvoDevoGeno/raw_reads_4_tony/RNA_trim"

# First, preview what would be generated
print("PREVIEW: Scanning for paired files...")
preview_count = preview_batch_generation(folder_path)

if preview_count > 0:
    response = input(f"\nProceed to generate {preview_count} batch files? (y/n): ").lower().strip()
    if response == 'y' or response == 'yes':
        print("\nGenerating batch files...")
        generate_rna_batch_files_and_master(folder_path, batch_output_folder, trim_output_folder)
    else:
        print("Batch file generation cancelled.")
else:
    print("No paired files found. Please check your file naming and folder path.")