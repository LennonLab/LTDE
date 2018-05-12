#!/bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=8,vmem=100gb,walltime=24:00:00
#PBS -M wrshoema@indiana.edu
#PBS -m abe
#PBS -j oe

# iRep and bPTR uses paired-end information, which breseq doesn't use
# so we'll need to map the reads using bwa and then clean them a little
# with samtools

module load bwa
module load samtools

mkdir -p /N/dc2/projects/muri2/Task2/LTDE/data/bwa_sam

for folder in "/N/dc2/projects/muri2/Task2/LTDE/data/reads_clean/"*
do
  R1="${folder}/"*_R1_paired.fastq.gz
  R2="${folder}/"*_R2_paired.fastq.gz
  strain="$(  echo "$folder" | cut -d"/" -f10-10 | cut -d"-" -f2-2 )"
  strain_rep="$(  echo "$folder" | cut -d"/" -f10-10 | cut -d"-" -f2-3 )"
  ref="/N/dc2/projects/muri2/Task2/LTDE/data/reference_genomes/"*"/${strain}/G-Chr1.fna"
  #bwa index $ref
  #samtools faidx $ref
  #bwa_sam_out=/N/dc2/projects/muri2/Task2/LTDE/data/bwa_sam/${strain_rep}.sam
  #bwa mem -t 4 $ref $R1 $R2 | samtools view -F 4 -bT $ref - \
      | samtools sort -o - | samtools view -h -o $bwa_sam_out
done


ATCC43928_C="/N/dc2/projects/muri2/Task2/LTDE/data/bwa_sam/ATCC43928-C"*".sam"
ATCC43928_C_out="/N/dc2/projects/muri2/Task2/LTDE/data/bwa_sam/ATCC43928-C.sam"
samtools merge - $ATCC43928_C | samtools sort -o $ATCC43928_C_out
