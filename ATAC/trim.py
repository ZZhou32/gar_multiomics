import os

def get_file_bases(folder_path):
    """
    Extract the base of the file name from all .fastaq.gz files in the specified folder.
    The base is defined as everything before 'R1' or 'R2'.
    """
    file_bases = set()  # Use a set to avoid duplicates
    for file_name in os.listdir(folder_path):
        if file_name.endswith(".fastq.gz"):
            if "R1" in file_name:
                base = file_name.split("_R1")[0]
            elif "R2" in file_name:
                base = file_name.split("_R2")[0]
            else:
                continue  # Skip files that don't match the pattern
            file_bases.add(base)
    return file_bases


def generate_batch_files_and_master(folder_path, batch_folder,trim_folder):
    """
    Generate individual batch files for each base name and a master batch file to submit them.
    """
    file_bases = get_file_bases(folder_path)

    # Ensure the output folder exists
    os.makedirs(batch_folder, exist_ok=True)

    # Open the master batch file for writing
    master_batch_path = os.path.join(batch_folder, "master_trim.sb")
    with open(master_batch_path, "w") as master_batch_file:
        master_batch_file.write("#!/bin/bash --login\n\n")  # Header for the master batch file

        for base in file_bases:
            batch_file_name = os.path.join(batch_folder, f"{base}.sb")  # For Linux/MacOS

            # Customize the content of the batch file
            batch_content = f"""#!/bin/bash --login
            
echo "Processing prefix: {base}"
# Number of nodes
#SBATCH --nodes=8

# Number of tasks to run on each node
#SBATCH --ntasks-per-node=6

# Number of CPUs per task
#SBATCH --cpus-per-task=4

# Memory per Node
# Specify "M" or "G" for MB and GB respectively
#SBATCH --mem=8G

# Wall time
# Format: "minutes", "hours:minutes:seconds", 
# "days-hours", or "days-hours:minutes"
#SBATCH --time=03:00:00

# Mail address
#SBATCH --mail-user=zhouzeh2@msu.edu

# Standard output and error to file
# %x: job name, %j: job ID
#SBATCH --output=/mnt/research/FishEvoDevoGeno/raw_reads_4_tony/ATAC_trim/out/%x-%j.SLURMout

# Set the input directory dynamically based on the base name
INPUT_DIR="{folder_path}"
OUTPUT_DIR="{trim_folder}"
in1="${{INPUT_DIR}}/{base}_R1.fastq.gz"
in2="${{INPUT_DIR}}/{base}_R2.fastq.gz"
out1="${{OUTPUT_DIR}}/{base}_R1NOADAPTER.fq"
out2="${{OUTPUT_DIR}}/{base}_R2NOADAPTER.fq"

module purge ## to remove current modules
#module load iccifort/2020.1.217
module load  BBMap/39.01-GCC-12.3.0  

# Run bbduk.sh (ensure the command is correct for your use case)
bash bbduk.sh \\
    in1="$in1" \\
    in2="$in2" \\
    out1="$out1" \\
    out2="$out2" \\
    ref="/mnt/research/FishEvoDevoGeno/raw_reads_4_tony/adapter_contamination_sequences.fa" \\
    ktrim=rl k=23 mink=11 hdist=0 tpe tbo
"""
            # Write the batch file
            with open(batch_file_name, "w") as batch_file:
                batch_file.write(batch_content)

            # Add sbatch command to the master batch file
            master_batch_file.write(f"sbatch script/ATAC/trim/{os.path.basename(batch_file_name)}\n")

    print(f"Generated {len(file_bases)} batch files and a master batch file in {batch_folder}")


# Example usage
folder_path = "/mnt/ufs18/rs-032/FishEvoDevoGeno/raw_reads_4_tony/ATAC" 
batch_output_folder = "/mnt/ufs18/rs-032/FishEvoDevoGeno/raw_reads_4_tony/script/ATAC/trim" 
trim_output_folder = "/mnt/ufs18/rs-032/FishEvoDevoGeno/raw_reads_4_tony/ATAC_trim"

generate_batch_files_and_master(folder_path, batch_output_folder,trim_output_folder)

