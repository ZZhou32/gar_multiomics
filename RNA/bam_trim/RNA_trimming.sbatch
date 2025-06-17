#!/bin/bash --login
#SBATCH --job-name=st_21_22_1__bam_trim
#SBATCH --output=/mnt/research/FishEvoDevoGeno/raw_reads_4_tony/RNA_bam/out/st_21_22_1__%j.out
#SBATCH --error=/mnt/research/FishEvoDevoGeno/raw_reads_4_tony/RNA_bam/out/st_21_22_1__%j.err
#SBATCH --time=16:00:00
#SBATCH --mem=64GB
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10

module purge
module load SAMtools/1.18-GCC-12.3.0
module load picard/2.25.1-Java-11
module load parallel/20230722-GCCcore-12.3.0
module load deepTools/3.5.5-foss-2023a

cd /mnt/research/FishEvoDevoGeno/raw_reads_4_tony/RNA_bam

echo "Processing samples"
# Mark duplicates using Picard
parallel --jobs 3 'java -Xmx32G -jar $EBROOTPICARD/picard.jar MarkDuplicates -I {} -O {.}_no_dup.bam -M {.}_dups.txt -REMOVE_DUPLICATES true' ::: *__sorted.bam

# remove non-uniquely-mapped reads from all ATAC reads
parallel 'samtools view -b -q 10 {} > {.}_uniq.bam' ::: *__sorted_no_dup.bam
echo "Finished sorting"


# index filtered bam files for ATAC reads
parallel 'samtools index {}' ::: *__sorted_no_dup_uniq.bam
echo "Finished indexing"