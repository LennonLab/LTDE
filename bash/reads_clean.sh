#!/bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=8,vmem=100gb,walltime=24:00:00
#PBS -M wrshoema@indiana.edu
#PBS -m abe
#PBS -j oe

module load fastqc
module load cutadapt

# trim data and remove adaptors.
mkdir -p /N/dc2/projects/muri2/Task2/LTDE/data/reads_clean
mkdir -p /N/dc2/projects/muri2/Task2/LTDE/data/reads_clean_quality

for R1 in "/N/dc2/projects/muri2/Task2/LTDE/data/reads_raw/"*"_R1_"*".fastq.gz"
do
  R2="${R1/_R1_/_R2_}"
  R1_clean="${R1/.fastq.gz/_clean.fastq.gz}"
  R1_clean="${R1_clean/reads_raw/reads_clean}"
  R2_clean="${R2/.fastq.gz/_clean.fastq.gz}"
  R2_clean="${R2_clean/reads_raw/reads_clean}"

  cutadapt -q 15,15 -u 10 -u -5 \
            -a file:/N/dc2/projects/muri2/Task2/LTDE/bash/transposase.fa \
            --minimum-length 20 \
            -o $R1_clean \
            -p $R2_clean \
            $R1 $R2

  fastqc_R1="${R1_clean/reads_clean/reads_clean_quality}"
  fastqc_R1="$( echo "${fastqc_R1}" | cut -d "." -f1-1 )"
  fastqc_R2="${R2_clean/reads_clean/reads_clean_quality}"
  fastqc_R2="$( echo "${fastqc_R2}" | cut -d "." -f1-1 )"
  mkdir $fastqc_R1
  mkdir $fastqc_R2
  fastqc $R1_clean --outdir=$fastqc_R1
  fastqc $R2_clean --outdir=$fastqc_R2
done
