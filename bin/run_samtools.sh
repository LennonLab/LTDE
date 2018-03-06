#!/bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=8,vmem=100gb,walltime=48:00:00
#PBS -M wrshoema@indiana.edu
#PBS -m abe
#PBS -j oe

module load bwa
module load samtools

out_sam="/N/dc2/projects/muri2/Task2/LTDE/data/sam_output/"
mkdir -p $out_sam

for folder in "/N/dc2/projects/muri2/Task2/LTDE/data/reads_clean/"*
do
  R1="${folder}/"*_R1_paired.fastq.gz
  R2="${folder}/"*_R2_paired.fastq.gz
  strain="$(  echo "$folder" | cut -d"/" -f10-10 | cut -d"-" -f2-2 )"
  strain_rep="$(  echo "$folder" | cut -d"/" -f10-10 | cut -d"-" -f2-3 )"
  #if [ "${strain_rep}" = "ATCC43928-C1" ]; then
  #  C1="/N/dc2/projects/muri2/Task2/LTDE/data/reads_clean/GSF1046-ATCC43928-C2/GSF1046-ATCC43928-C2_clean_R1_paired.fastq.gz"
  #  C2="/N/dc2/projects/muri2/Task2/LTDE/data/reads_clean/GSF1046-ATCC43928-C2/GSF1046-ATCC43928-C2_clean_R2_paired.fastq.gz"
  #  reads="${R1} ${R2} ${C1} ${C2}"
  #  out_sam_strain="${out_sam}ATCC43928-C"
  #elif [ "${strain_rep}" = "ATCC43928-C2" ]; then
  #  continue
  #else
  #  continue
  out_sam_strain="${out_sam}${strain_rep}"
  out_sam_strain="${out_sam}${strain_rep}"
  #fi
  ref="/N/dc2/projects/muri2/Task2/LTDE/data/reference_genomes/"*"/${strain}/G-Chr1.fna"
  bwa index $ref
  samtools faidx $ref
  mkdir -p $out_sam_strain
  mkdir -p "${out_sam_strain}/tmp"
  #bwa mem -t 4 $ref $R1 $R2 > "${out_sam_strain}/${strain_rep}.sam"
  # mapped reads
  #samtools view -F 4 -bT $ref "${out_sam_strain}/${strain_rep}.sam" > "${out_sam_strain}/${strain_rep}_mapped.bam"
  #samtools sort -T "${out_sam_strain}/tmp" -o "${out_sam_strain}/${strain_rep}_mapped_sort.bam" "${out_sam_strain}/${strain_rep}_mapped.bam"
  #samtools index "${out_sam_strain}/${strain_rep}_mapped_sort.bam"
  #samtools rmdup "${out_sam_strain}/${strain_rep}_mapped_sort.bam" "${out_sam_strain}/${strain_rep}_mapped_sort_NOdup.bam"
  #samtools index "${out_sam_strain}/${strain_rep}_mapped_sort_NOdup.bam"
  samtools sort -T "${out_sam_strain}/tmp" -o "${out_sam_strain}/${strain_rep}_mapped_sort_NOdup_sort.bam" "${out_sam_strain}/${strain_rep}_mapped_sort_NOdup.bam" 
  samtools index "${out_sam_strain}/${strain_rep}_mapped_sort_NOdup_sort.bam"
  samtools view -h -o "${out_sam_strain}/${strain_rep}_mapped_sort_NOdup_sort.sam" "${out_sam_strain}/${strain_rep}_mapped_sort_NOdup_sort.bam"
done


# then merge the C lines........
