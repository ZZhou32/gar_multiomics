#!/bin/bash

if [ $# -lt 4 ]
then
    echo ""
    echo "Usage: $0 <gar|hg38|hg19|mm10> <peak.file.lst> <bam.file.lst> <output prefix> [PEAKS|FOOTPRINTS]"
    echo ""
    exit -1
fi

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")

# Activate environment
source /mnt/research/FishEvoDevoGeno/raw_reads_4_tony/script/ATAC_analysis/counts.sh # Or your conda.sh path
conda activate atac

# CMD parameters
ATYPE=${1}
PFLST=${2}
BFLST=${3}
OP=${4}
FOOTPRINTS=0
if [ $# -eq 5 ]
then
    if [ ${5} == "FOOTPRINTS" ]
    then
        FOOTPRINTS=1
    fi
fi

# Gar-specific settings
if [ ${ATYPE} == "gar" ]
then
    PROMOTER_BED="${BASEDIR}/../bed/gar.promoter.bed"  # Uncompressed
    BLACKLIST="${BASEDIR}/../bed/gar_blacklist.bed"    # Create if available
    CHR_FILTER="grep -v \"^Un\" | grep -v \"^scaffold\""  # Filter unusual contigs
else
    PROMOTER_BED="${BASEDIR}/../bed/${ATYPE}.promoter.bed.gz"
    BLACKLIST="${BASEDIR}/../bed/wgEncodeDacMapabilityConsensusExcludable.bed.gz"
    CHR_FILTER="grep -v \"^Y\" | grep -v \"chrY\""
fi

# Concatenate peaks
rm -f ${OP}.concat.peaks
for P in `cat ${PFLST}`
do
    if [ ! -f ${P} ]
    then
        echo "Peak file does not exist:" ${P}
        exit -1
    else
        echo "Reading" ${P} ", Avg. peak width" `cat ${P} | awk '{SUM+=($3-$2);} END {print SUM/NR;}'`
        cat ${P} >> ${OP}.concat.peaks
    fi
done

# Cluster peaks
bedtools merge -i <(sort -k1,1V -k2,2n ${OP}.concat.peaks) | awk '{print $0"\tPeak"sprintf("%08d", NR);}' | gzip -c > ${OP}.clustered.peaks.gz
rm ${OP}.concat.peaks

# Remove blacklisted regions (skip if no gar blacklist exists)
if [ -f "${BLACKLIST}" ] || [ -f "${BLACKLIST}.gz" ]
then
    echo "Filtering blacklisted regions"
    if [[ ${BLACKLIST} == *.gz ]]
    then
        bedtools intersect -v -a <(zcat ${OP}.clustered.peaks.gz) -b <(zcat ${BLACKLIST}) | eval ${CHR_FILTER} | gzip -c > ${OP}.clustered.peaks.gz.tmp
    else
        bedtools intersect -v -a <(zcat ${OP}.clustered.peaks.gz) -b <(cat ${BLACKLIST}) | eval ${CHR_FILTER} | gzip -c > ${OP}.clustered.peaks.gz.tmp
    fi
    mv ${OP}.clustered.peaks.gz.tmp ${OP}.clustered.peaks.gz
else
    echo "No blacklist found for ${ATYPE}, skipping filtering"
    zcat ${OP}.clustered.peaks.gz | eval ${CHR_FILTER} | gzip -c > ${OP}.clustered.peaks.gz.tmp
    mv ${OP}.clustered.peaks.gz.tmp ${OP}.clustered.peaks.gz
fi

# Check post-clustered peak width
echo "Post-clustered peak width" `zcat ${OP}.clustered.peaks.gz | awk '{SUM+=($3-$2);} END {print SUM/NR;}'`

# Count fragments in peaks
for B in `cat ${BFLST}`
do
    if [ ! -f ${B} ]
    then
        echo "BAM file does not exist:" ${B}
        exit -1
    else
        if [ ${FOOTPRINTS} -eq 1 ]
        then
            alfred count_dna -f 200,1000 -o ${OP}.count.gz -i ${OP}.clustered.peaks.gz ${B}
        else
            alfred count_dna -o ${OP}.count.gz -i ${OP}.clustered.peaks.gz ${B}
        fi
        SID=`zcat ${OP}.count.gz | head -n 1 | cut -f 5`
        mv ${OP}.count.gz ${OP}.${SID}.count.gz
    fi
done
rm ${OP}.clustered.peaks.gz

# Create count matrix
rm -f ${OP}.counts
for F in ${OP}.*.count.gz
do
    echo "Aggregating" ${F}
    if [ ! -f ${OP}.counts ]
    then
        zcat ${F} > ${OP}.counts
    else
        paste ${OP}.counts <(zcat ${F} | cut -f 5) > ${OP}.counts.tmp
        mv ${OP}.counts.tmp ${OP}.counts
    fi
    rm ${F}
done
gzip ${OP}.counts

# Split counts into TSS peaks and non-TSS peaks
if [ -f "${PROMOTER_BED}" ] || [ -f "${PROMOTER_BED}.gz" ]
then
    echo "Splitting TSS/non-TSS peaks"
    zcat ${OP}.counts.gz | head -n 1 > ${OP}.tss.counts
    if [[ ${PROMOTER_BED} == *.gz ]]
    then
        bedtools intersect -wa -a <(zcat ${OP}.counts.gz | tail -n +2) -b <(zcat ${PROMOTER_BED}) | sort -k1,1V -k2,2n | uniq >> ${OP}.tss.counts
    else
        bedtools intersect -wa -a <(zcat ${OP}.counts.gz | tail -n +2) -b <(cat ${PROMOTER_BED}) | sort -k1,1V -k2,2n | uniq >> ${OP}.tss.counts
    fi
    
    zcat ${OP}.counts.gz | head -n 1 > ${OP}.notss.counts
    if [[ ${PROMOTER_BED} == *.gz ]]
    then
        bedtools intersect -v -wa -a <(zcat ${OP}.counts.gz | tail -n +2) -b <(zcat ${PROMOTER_BED}) | sort -k1,1V -k2,2n | uniq >> ${OP}.notss.counts
    else
        bedtools intersect -v -wa -a <(zcat ${OP}.counts.gz | tail -n +2) -b <(cat ${PROMOTER_BED}) | sort -k1,1V -k2,2n | uniq >> ${OP}.notss.counts
    fi
    gzip ${OP}.tss.counts ${OP}.notss.counts
else
    echo "No promoter file found for ${ATYPE}, skipping TSS splitting"
fi

# Deactivate environment
source deactivate