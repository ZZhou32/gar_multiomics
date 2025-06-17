#!/bin/bash --login
#SBATCH --job-name=st_22_23_3__align
#SBATCH --output=/mnt/research/FishEvoDevoGeno/raw_reads_4_tony/script/RNA/mapping/st_22_23_3__%j.out
#SBATCH --error=/mnt/research/FishEvoDevoGeno/raw_reads_4_tony/script/RNA/mapping/st_22_23_3__%j.err
#SBATCH --time=04:00:00
#SBATCH --mem=16GB
#SBATCH --cpus-per-task=4

module purge

module load Bowtie2/2.5.1-GCC-12.3.0
module load SAMtools/1.18-GCC-12.3.0

echo "Processing sample: st_22_23_3_"
bowtie2 --very-sensitive -x /mnt/research/FishEvoDevoGeno/raw_reads_4_tony/index_genome/fLepOcu1hap2 -1 /mnt/research/FishEvoDevoGeno/raw_reads_4_tony/RNA_trim/st_22_23_3__1NOADAPTERR.fq.gz -2 /mnt/research/FishEvoDevoGeno/raw_reads_4_tony/RNA_trim/st_22_23_3__2NOADAPTERR.fq.gz | samtools view -u | samtools sort -o /mnt/research/FishEvoDevoGeno/raw_reads_4_tony/RNA_bam/st_22_23_3__sorted.bam 
echo "Finished alignment for st_22_23_3_"
