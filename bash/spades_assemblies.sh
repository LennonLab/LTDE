#!/bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=8,vmem=100gb,walltime=01:00:00
#PBS -M wrshoema@iu.edu
#PBS -m abe
#PBS -j oe

module load spades

ltde=/N/dc2/projects/muri2/Task2/LTDE





# run 3

#spades.py --careful --cov-cutoff auto \
#      --nanopore "${ltde}/data/nanopore_basecalled_bc_merged/ATCC13985_clean.fastq" \
#      --pe1-1 "${ltde}/illumina_data/ATCC13985_D400_100/Ev13483_TAAGGCGAATCT-AAGGAGTA_L002_R1_001_clean.fastq.gz" \
#      --pe1-2 "${ltde}/illumina_data/ATCC13985_D400_100/Ev13483_TAAGGCGAATCT-AAGGAGTA_L002_R2_001_clean.fastq.gz" \
#      -o "${ltde}/data/spades_assemblies/ATCC13985"

spades.py --careful --cov-cutoff auto \
      --nanopore "${ltde}/data/nanopore_basecalled_bc_merged/ATCC43928_clean.fastq" \
      --pe1-1 "${ltde}/illumina_data/ATCC43928_D400_100/Ev43828_TAAGGCGAATCT-ACTGCATA_L002_R1_001_clean.fastq.gz" \
      --pe1-2 "${ltde}/illumina_data/ATCC43928_D400_100/Ev43828_TAAGGCGAATCT-ACTGCATA_L002_R2_001_clean.fastq.gz" \
      -o "${ltde}/data/spades_assemblies/ATCC43928"

#spades.py --careful --cov-cutoff auto \
#      --nanopore "${ltde}/data/nanopore_basecalled_bc_merged/KBS0701_clean.fastq" \
#      --pe1-1 "${ltde}/illumina_data/KBS0701_2015_SoilGenomes/KBS0701_R1_clean.fastq" \
#      --pe1-2 "${ltde}/illumina_data/KBS0701_2015_SoilGenomes/KBS0701_R2_clean.fastq" \
#      -o "${ltde}/data/spades_assemblies/KBS0701"

spades.py --careful --cov-cutoff auto \
      --nanopore "${ltde}/data/nanopore_basecalled_bc_merged/KBS0702_clean.fastq" \
      --pe1-1 "${ltde}/illumina_data/KBS0702_2015_SoilGenomes/KBS0702_R1_clean.fastq" \
      --pe1-2 "${ltde}/illumina_data/KBS0702_2015_SoilGenomes/KBS0702_R2_clean.fastq" \
      --pe2-1 "${ltde}/illumina_data/KBS0702_D400_100/Ev702_TAAGGCGAATCT-TAGATCGC_L002_R1_001_clean.fastq.gz" \
      --pe2-2 "${ltde}/illumina_data/KBS0702_D400_100/Ev702_TAAGGCGAATCT-TAGATCGC_L002_R2_001_clean.fastq.gz" \
      -o "${ltde}/data/spades_assemblies/KBS0702"

spades.py --careful --cov-cutoff auto \
      --nanopore "${ltde}/data/nanopore_basecalled_bc_merged/KBS0703_clean.fastq" \
      --pe1-1 "${ltde}/illumina_data/KBS0703_GSF911/GSF911-Ar_S10_L001_R1_001_clean.fastq.gz" \
      --pe1-2 "${ltde}/illumina_data/KBS0703_GSF911/GSF911-Ar_S10_L001_R2_001_clean.fastq.gz" \
      -o "${ltde}/data/spades_assemblies/KBS0703"

spades.py --careful --cov-cutoff auto \
      --nanopore "${ltde}/data/nanopore_basecalled_bc_merged/KBS0705_clean.fastq" \
      --pe1-1 "${ltde}/illumina_data/KBS0705_2015_SoilGenomes/KBS0705_R1_clean.fastq" \
      --pe1-2 "${ltde}/illumina_data/KBS0705_2015_SoilGenomes/KBS0705_R2_clean.fastq" \
      -o "${ltde}/data/spades_assemblies/KBS0705"

spades.py --careful --cov-cutoff auto \
      --nanopore "${ltde}/data/nanopore_basecalled_bc_merged/KBS0707_clean.fastq" \
      --pe1-1 "${ltde}/illumina_data/KBS0707_D400_100/Ev707_TAAGGCGAATCT-CTCTCTAT_L002_R1_001_clean.fastq.gz" \
      --pe1-2 "${ltde}/illumina_data/KBS0707_D400_100/Ev707_TAAGGCGAATCT-CTCTCTAT_L002_R2_001_clean.fastq.gz" \
      -o "${ltde}/data/spades_assemblies/KBS0707"

spades.py --careful --cov-cutoff auto \
      --nanopore "${ltde}/data/nanopore_basecalled_bc_merged/KBS0724_clean.fastq" \
      --pe1-1 "${ltde}/illumina_data/KBS0724_2015_SoilGenomes/KBS0724_R1_clean.fastq" \
      --pe1-2 "${ltde}/illumina_data/KBS0724_2015_SoilGenomes/KBS0724_R2_clean.fastq" \
      -o "${ltde}/data/spades_assemblies/KBS0724"


#prokka --compliant --centre IUB --outdir PRJNA539822 --locustag FCE86 \
#        --prefix FCE86-Genome \
#        contigs.fa
