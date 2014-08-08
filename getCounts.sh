#!/bin/bash
#$ -cwd

bam=$1
annotation=$2
sample=$3

tempCount=$bam.tempCount
countFile=$bam.counts.txt
# samtools view -h $bam | htseq-count -s no - $annotation >> $tempCount

if [[ $sample == "" ]]; then
    pwd=`pwd`
    sample=`basename $pwd | cut -f1 -d '.'`
fi

featureCounts -t exon -g gene_id -a $annotation -o $tempCount $bam

echo -e "gene\t$sample" > $countFile

cat $tempCount | awk 'NR > 2 {print $1"\t"$7}' >> $countFile

# rm $tempCount