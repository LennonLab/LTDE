#!/bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=8,vmem=50gb,walltime=2:00:00
#PBS -M wrshoema@indiana.edu
#PBS -m abe
#PBS -j oe

module load bwa
module load samtools
module load breseq
#for folder in "/N/dc2/projects/muri2/Task2/LTDE/data/reads_clean/"*
#do
#  R1="${folder}/"*_R1_paired.fastq.gz
#  R2="${folder}/"*_R2_paired.fastq.gz

R1="/N/dc2/projects/muri2/Task2/LTDE/data/reads_clean/GSF1046-ATCC13985-A/GSF1046-ATCC13985-A_clean_R1_paired.fastq.gz"
R2="/N/dc2/projects/muri2/Task2/LTDE/data/reads_clean/GSF1046-ATCC13985-A/GSF1046-ATCC13985-A_clean_R2_paired.fastq.gz"
#ref="/N/dc2/projects/muri2/Task2/LTDE/data/reference_genomes/2016_KBSGenomes_Annotate/ATCC13985/G-Chr1.fna"
ref='/N/dc2/projects/muri2/Task2/LTDE/data/reference_genomes/genomes_rename_fna/ATCC13985.fna'
#bwa index $ref
#samtools faidx $ref

out_rmdup_test="/N/dc2/projects/muri2/Task2/LTDE/data/bwa_output/ATCC13985-A/test_rmdup.bam"
out_rmdup_sort_test="/N/dc2/projects/muri2/Task2/LTDE/data/bwa_output/ATCC13985-A/test_rmdup_sort.bam"
out_sam="/N/dc2/projects/muri2/Task2/LTDE/data/bwa_output/ATCC13985-A/out_sam.sam"

out_noredup_sam="/N/dc2/projects/muri2/Task2/LTDE/data/bwa_output/ATCC13985-A/out_noredup_sam.sam"
#bwa mem -t 4 $ref $R1 $R2 | samtools view -F 4 -bT $ref | samtools sort -o $out_test
#bwa mem -t 4 $ref $R1 $R2 | samtools view -F 4 -bT $ref - | samtools sort -o - | samtools rmdup - - | samtools view -h -o $out_sam
#bwa mem -t 4 $ref $R1 $R2 | samtools view -F 4 -bT $ref - | samtools sort -o - | samtools rmdup - -| samtools sort -o $out_rmdup_sort_test

#bwa mem -t 4 $ref $R1 $R2 | samtools view -F 4 -bT $ref - \
#    | samtools sort -o - | samtools rmdup - -| samtools sort -n -o - \
#    | samtools view -h -o $out_sam

# do GATK alignment?

bwa mem -t 4 $ref $R1 $R2 | samtools view -F 4 -bT $ref - \
    | samtools sort -o - | samtools view -h -o $out_noredup_sam

OUT_breseq_strain="/N/dc2/projects/muri2/Task2/LTDE/data/bwa_output/ATCC13985-A/breseq_noredup_test"
breseq -j 8 -p -o $OUT_breseq_strain -r $ref --aligned-sam $out_noredup_sam
