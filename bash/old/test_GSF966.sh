#!/bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=8,vmem=100gb,walltime=18:00:00
#PBS -M wrshoema@iu.edu
#PBS -m abe
#PBS -j oe

#module load cutadapt
#module load spades
module load unicycler

GSF_PATH=/N/dc2/projects/muri2/Task2/LTDE/GSF966

#cutadapt -q 30,30 -u 20 \
#          -b AGATCGGAAGAGCACACGTCTGAACTCCAGTCACACAGTGATCTCGTATG \
#          -b AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCC \
#          --pair-filter=any \
#          --minimum-length 20 \
#          -n 2 \
#          -o "${GSF_PATH}/GSF966-1-Arthro-6k_S1_R1_001_clean.fastq" \
#          -p "${GSF_PATH}/GSF966-1-Arthro-6k_S1_R2_001_clean.fastq" \
#          "${GSF_PATH}/GSF966-1-Arthro-6k_S1_R1_001.fastq.gz" \
#          "${GSF_PATH}/GSF966-1-Arthro-6k_S1_R2_001.fastq.gz"


#cutadapt -q 30,30 -u 20 \
#          -b AGATCGGAAGAGCACACGTCTGAACTCCAGTCACACAGTGATCTCGTATG \
#          -b AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCC \
#          --pair-filter=any \
#          --minimum-length 20 \
#          -n 2 \
#          -o "${GSF_PATH}/GSF966-2-Arthro-13k_S2_R1_001_clean.fastq" \
#          -p "${GSF_PATH}/GSF966-2-Arthro-13k_S2_R2_001_clean.fastq" \
#          "${GSF_PATH}/GSF966-2-Arthro-13k_S2_R1_001.fastq.gz" \
#          "${GSF_PATH}/GSF966-2-Arthro-13k_S2_R2_001.fastq.gz"



#spades.py --careful --cov-cutoff auto \
#      --nanopore /N/dc2/projects/muri2/Task2/LTDE/data/nanopore_basecalled_bc_merged/KBS0703_clean.fastq \
#      --pe1-1 "${GSF_PATH}/GSF966-1-Arthro-6k_S1_R1_001_clean.fastq" \
#      --pe1-2 "${GSF_PATH}/GSF966-1-Arthro-6k_S1_R2_001_clean.fastq" \
#      --pe2-1 "${GSF_PATH}/GSF966-2-Arthro-13k_S2_R1_001_clean.fastq" \
#      --pe2-2 "${GSF_PATH}/GSF966-2-Arthro-13k_S2_R2_001_clean.fastq" \
#      -o /N/dc2/projects/muri2/Task2/LTDE/GSF966_KBS0703_spades_13k_6k


#spades.py --careful --cov-cutoff auto \
#      --nanopore /N/dc2/projects/muri2/Task2/LTDE/data/nanopore_basecalled_bc_merged/KBS0703_clean.fastq \
#      --pe1-1 "${GSF_PATH}/GSF966-1-Arthro-6k_S1_R1_001_clean.fastq" \
#      --pe1-2 "${GSF_PATH}/GSF966-1-Arthro-6k_S1_R2_001_clean.fastq" \
#      -o /N/dc2/projects/muri2/Task2/LTDE/GSF966_KBS0703_spades_6k


#spades.py --careful --cov-cutoff auto \
#      --nanopore /N/dc2/projects/muri2/Task2/LTDE/data/nanopore_basecalled_bc_merged/KBS0703_clean.fastq \
#      --pe1-1 "${GSF_PATH}/GSF966-2-Arthro-13k_S2_R1_001_clean.fastq" \
#      --pe1-2 "${GSF_PATH}/GSF966-2-Arthro-13k_S2_R2_001_clean.fastq" \
#      -o /N/dc2/projects/muri2/Task2/LTDE/GSF966_KBS0703_spades_13k



unicycler \
    -1 "${GSF_PATH}/GSF966-1-Arthro-6k_S1_R1_001_clean.fastq" \
    -2 "${GSF_PATH}/GSF966-1-Arthro-6k_S1_R2_001_clean.fastq" \
    -l /N/dc2/projects/muri2/Task2/LTDE/data/nanopore_basecalled_bc_merged/KBS0703_clean.fastq \
    -o /N/dc2/projects/muri2/Task2/LTDE/GSF966_KBS0703_uni_6k


unicycler \
    -1 "${GSF_PATH}/GSF966-2-Arthro-13k_S2_R1_001_clean.fastq" \
    -2 "${GSF_PATH}/GSF966-2-Arthro-13k_S2_R2_001_clean.fastq" \
    -l /N/dc2/projects/muri2/Task2/LTDE/data/nanopore_basecalled_bc_merged/KBS0703_clean.fastq \
    -o /N/dc2/projects/muri2/Task2/LTDE/GSF966_KBS0703_uni_13k
