import os

def load_base_names_from_file(base_names_file_path, folder_path):
    """
    Load base names from RNA_base_names.txt file and verify corresponding files exist.
    """
    if not os.path.exists(base_names_file_path):
        print(f"Error: RNA_base_names.txt file not found at {base_names_file_path}")
        return {}
    
    if not os.path.exists(folder_path):
        print(f"Error: Folder {folder_path} does not exist!")
        return {}
    
    # Read base names from file
    with open(base_names_file_path, 'r') as f:
        base_names = [line.strip() for line in f if line.strip()]
    
    print(f"Loaded {len(base_names)} base names from RNA_base_names.txt")
    
    # Dictionary to store file pairs: {base_name: {'R1': file1, 'R2': file2}}
    paired_files = {}
    
    # For each base name, verify the corresponding files exist
    for base_name in base_names:
        print(f"Processing: {base_name}")
        
        # Construct expected file names (assuming trimmed files with NOADAPTERR suffix)
        r1_file = f"{base_name}_1NOADAPTERR.fq.gz"
        r2_file = f"{base_name}_2NOADAPTERR.fq.gz"
        
        r1_path = os.path.join(folder_path, r1_file)
        r2_path = os.path.join(folder_path, r2_file)
        
        # Check if both files exist
        if os.path.exists(r1_path) and os.path.exists(r2_path):
            paired_files[base_name] = {
                'R1': r1_file,
                'R2': r2_file,
                'R1_path': r1_path,
                'R2_path': r2_path
            }
            print(f"  ✓ Found pair: {base_name}")
            print(f"    R1: {r1_file}")
            print(f"    R2: {r2_file}")
        else:
            print(f"  ✗ Missing files for {base_name}")
            if not os.path.exists(r1_path):
                print(f"    Missing R1: {r1_file}")
            if not os.path.exists(r2_path):
                print(f"    Missing R2: {r2_file}")
    
    print(f"\nSummary:")
    print(f"Base names from file: {len(base_names)}")
    print(f"Complete pairs found: {len(paired_files)}")
    print(f"Missing pairs: {len(base_names) - len(paired_files)}")
    
    if paired_files:
        print(f"\nPaired samples found:")
        for base_name in sorted(paired_files.keys()):
            print(f"  - {base_name}")
    
    return paired_files

def generate_salmon_batch_files_and_master(folder_path, batch_folder, salmon_index, quant_folder, base_names_file_path):
    """
    Generate individual batch files for Salmon quantification based on base names from RNA_base_names.txt file.
    """
    # Load base names from file and verify files exist
    paired_files = load_base_names_from_file(base_names_file_path, folder_path)
    
    if not paired_files:
        print("No paired files found! Cannot generate batch scripts.")
        return

    # Ensure the output folders exist
    os.makedirs(batch_folder, exist_ok=True)
    os.makedirs(quant_folder, exist_ok=True)
    os.makedirs(f"{quant_folder}/out", exist_ok=True)
    print(f"\nCreated output directories:")
    print(f"  - Batch files: {batch_folder}")
    print(f"  - Quantification output: {quant_folder}")
    print(f"  - Logs: {quant_folder}/out")

    # Open the master batch file for writing
    master_batch_path = os.path.join(batch_folder, "master_salmon_quant.sb")
    
    with open(master_batch_path, "w") as master_batch_file:
        master_batch_file.write("#!/bin/bash --login\n\n")
        master_batch_file.write("# Master script to submit all Salmon quantification jobs\n")
        master_batch_file.write(f"# Generated for {len(paired_files)} paired samples\n")
        master_batch_file.write(f"# Source folder: {folder_path}\n\n")

        batch_files_created = 0
        
        for base_name in sorted(paired_files.keys()):
            file_info = paired_files[base_name]
            
            # Create batch file name using base name
            batch_file_name = os.path.join(batch_folder, f"{base_name}_salmon.sb")
            
            print(f"Creating batch file for: {base_name}")

            # Customize the content of the batch file for Salmon quantification - CORRECTED VERSION
            batch_content = f"""#!/bin/bash --login

#SBATCH --job-name={base_name}_salmon_job
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --mail-user=zhouzeh2@msu.edu
#SBATCH --output={quant_folder}/out/%x-%j.SLURMout
#SBATCH --error={quant_folder}/out/%x-%j.SLURMerr

echo "Processing: {base_name}"
echo "Job started at: $(date)"

# File paths
samp="{base_name}"
fn="{folder_path}"
in1="{file_info['R1_path']}"
in2="{file_info['R2_path']}"
salmon_index="{salmon_index}"
output_dir="{quant_folder}/{base_name}_quant"

# Create output directory
mkdir -p {quant_folder}
mkdir -p "$output_dir"

echo "Input file sizes:"
ls -lh "$in1" "$in2"

echo "Memory before processing:"
free -h

# Activate conda environment
conda activate salmon

echo "Processing sample ${{samp}}"
echo "Using Salmon version:"
salmon --version

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

if [[ ! -d "$salmon_index" ]]; then
    echo "ERROR: Salmon index not found: $salmon_index"
    exit 1
fi

echo "✓ Both input files and index verified"

# Run Salmon quantification
salmon quant -i "$salmon_index" -l A \\
    -1 "$in1" \\
    -2 "$in2" \\
    -p 8 --validateMappings -o "$output_dir"

if [[ $? -eq 0 ]]; then
    echo "✓ Salmon quantification completed successfully"
    echo "Output files:"
    ls -lh "$output_dir"
    
    if [[ -f "$output_dir/quant.sf" ]]; then
        echo "✓ Quantification file created"
        TRANSCRIPT_COUNT=$(tail -n +2 "$output_dir/quant.sf" | wc -l)
        echo "Transcripts quantified: $TRANSCRIPT_COUNT"
    fi
else
    echo "✗ Salmon quantification failed"
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
            master_batch_file.write(f"echo 'Submitting {base_name} Salmon quantification job...'\n")
            master_batch_file.write(f"sbatch {batch_file_name}\n")
            master_batch_file.write(f"echo 'Job submitted for {base_name}'\n")
            master_batch_file.write(f"sleep 2  # Delay between submissions\n\n")

        # Add summary to master batch file
        master_batch_file.write(f"echo 'All {batch_files_created} Salmon quantification jobs submitted!'\n")
        master_batch_file.write(f"echo 'Monitor with: squeue -u $USER'\n")
        master_batch_file.write(f"echo 'Check logs in: {quant_folder}/out/'\n")

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
    print(f"sbatch {batch_folder}/<sample_name>_salmon.sb")

def preview_batch_generation(folder_path, base_names_file_path):
    """
    Preview what batch files would be generated without actually creating them.
    """
    print("="*60)
    print("PREVIEW MODE - No files will be created")
    print("="*60)
    
    paired_files = load_base_names_from_file(base_names_file_path, folder_path)
    
    if paired_files:
        print(f"\nWould generate {len(paired_files)} batch files:")
        for base_name in sorted(paired_files.keys()):
            file_info = paired_files[base_name]
            print(f"  {base_name}_salmon.sb")
            print(f"    R1: {file_info['R1']}")
            print(f"    R2: {file_info['R2']}")
            print()
        
        # Write base names to a file during preview
        base_names_file = os.path.join(os.path.dirname(folder_path), "sample_base_names_preview.txt")
        with open(base_names_file, "w") as f:
            for base_name in sorted(paired_files.keys()):
                f.write(f"{base_name}\n")
        print(f"Base names written to: {base_names_file}")
    
    return len(paired_files)

# CORRECTED Configuration - Updated to match your working script
folder_path = "/mnt/research/FishEvoDevoGeno/raw_reads_4_tony/RNA_trim" 
batch_output_folder = "/mnt/research/FishEvoDevoGeno/raw_reads_4_tony/script/RNA/salmon" 
salmon_index = "/mnt/research/FishEvoDevoGeno/raw_reads_4_tony/salmon_prep/salmon_index"  # CORRECTED PATH
quant_output_folder = "/mnt/research/FishEvoDevoGeno/raw_reads_4_tony/salmon_quant"
base_names_file = "/mnt/research/FishEvoDevoGeno/raw_reads_4_tony/RNA_base_names.txt"

# First, preview what would be generated
print("PREVIEW: Loading base names from RNA_base_names.txt...")
preview_count = preview_batch_generation(folder_path, base_names_file)

if preview_count > 0:
    response = input(f"\nProceed to generate {preview_count} batch files? (y/n): ").lower().strip()
    if response == 'y' or response == 'yes':
        print("\nGenerating batch files...")
        generate_salmon_batch_files_and_master(folder_path, batch_output_folder, salmon_index, quant_output_folder, base_names_file)
    else:
        print("Batch file generation cancelled.")
else:
    print("No paired files found. Please check your RNA_base_names.txt file and folder path.")