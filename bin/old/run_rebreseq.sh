#!/bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=8,vmem=50gb,walltime=0:30:00
#PBS -M wrshoema@indiana.edu
#PBS -m abe
#PBS -j oe


# get uniquely mapped reads
# then re-do breseq alignment using sam file

#module load samtools
module load breseq

#breseq_bam="/N/dc2/projects/muri2/Task2/LTDE/data/breseq_output_v32/ATCC13985-A/data/reference.bam"
#breseq_sam="/N/dc2/projects/muri2/Task2/LTDE/data/rebreseq_output_v32/ATCC13985-A/reference.sam"
#breseq_header_sam="/N/dc2/projects/muri2/Task2/LTDE/data/rebreseq_output_v32/ATCC13985-A/reference_header.sam"
#breseq_sam_unique="/N/dc2/projects/muri2/Task2/LTDE/data/rebreseq_output_v32/ATCC13985-A/reference_unique.sam"

#samtools view -H $breseq_bam > $breseq_header_sam
#samtools view -F 4 $breseq_bam | grep "X1:i:1" | cat $breseq_header_sam - | \
#  samtools view -h -o $breseq_sam_unique

# now, re-run breseq
