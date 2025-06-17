#!/bin/bash --login
#SBATCH --job-name=submit_all
#SBATCH --output=/mnt/research/FishEvoDevoGeno/raw_reads_4_tony/script/RNA/mapping/submit_all_%j.out
#SBATCH --error=/mnt/research/FishEvoDevoGeno/raw_reads_4_tony/script/RNA/mapping/submit_all_%j.err
#SBATCH --time=5:00
#SBATCH --mem=1G
#SBATCH --partition=debug  # Adjust if needed

# Submit all individual jobs
for job in "/mnt/research/FishEvoDevoGeno/raw_reads_4_tony/script/RNA/mapping"/*RNA_mapping.sbatch; do
    if [[ "$job" != "/mnt/research/FishEvoDevoGeno/raw_reads_4_tony/script/RNA/mapping/submit_all.sbatch" ]]; then  # Avoid submitting itself
        echo "Submitting $job"
        sbatch "$job"
    fi
done
echo "All jobs submitted!"
