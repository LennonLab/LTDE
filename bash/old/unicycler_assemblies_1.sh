#!/bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=8,vmem=100gb,walltime=120:00:00
#PBS -M wrshoema@iu.edu
#PBS -m abe
#PBS -j oe

module load unicycler

ltde=/N/dc2/projects/muri2/Task2/LTDE



unicycler \
    -1 "${ltde}/illumina_data/ATCC13985_D400_100/Ev13483_TAAGGCGAATCT-AAGGAGTA_L002_R1_001_clean.fastq" \
    -2 "${ltde}/illumina_data/ATCC13985_D400_100/Ev13483_TAAGGCGAATCT-AAGGAGTA_L002_R2_001_clean.fastq" \
    -l "${ltde}/data/nanopore_basecalled_bc_merged/ATCC13985_clean.fastq" \
    -o "${ltde}/data/unicycler_assemblies/ATCC13985"


unicycler \
    -1 "${ltde}/illumina_data/ATCC43928_D400_100/Ev43828_TAAGGCGAATCT-ACTGCATA_L002_R1_001_clean.fastq" \
    -2 "${ltde}/illumina_data/ATCC43928_D400_100/Ev43828_TAAGGCGAATCT-ACTGCATA_L002_R1_001_clean.fastq" \
    -l "${ltde}/data/nanopore_basecalled_bc_merged/ATCC43928_clean.fastq" \
    -o "${ltde}/data/unicycler_assemblies/ATCC43928"


unicycler \
    -1 "${ltde}/illumina_data/KBS0701_2015_SoilGenomes/KBS0701_R1_clean.fastq" \
    -2 "${ltde}/illumina_data/KBS0701_2015_SoilGenomes/KBS0701_R2_clean.fastq" \
    -l "${ltde}/data/nanopore_basecalled_bc_merged/KBS0701_clean.fastq" \
    -o "${ltde}/data/unicycler_assemblies/KBS0701"


unicycler \
    -1 "${ltde}/illumina_data/KBS0702_2015_SoilGenomes/KBS0702_R1_clean.fastq" \
    -2 "${ltde}/illumina_data/KBS0702_2015_SoilGenomes/KBS0702_R2_clean.fastq" \
    -l "${ltde}/data/nanopore_basecalled_bc_merged/KBS0702_clean.fastq" \
    -o "${ltde}/data/unicycler_assemblies/KBS0702_1"


unicycler \
    -1 "${ltde}/illumina_data/KBS0702_D400_100/Ev702_TAAGGCGAATCT-TAGATCGC_L002_R1_001_clean.fastq.gz" \
    -2 "${ltde}/illumina_data/KBS0702_D400_100/Ev702_TAAGGCGAATCT-TAGATCGC_L002_R2_001_clean.fastq.gz" \
    -l "${ltde}/data/nanopore_basecalled_bc_merged/KBS0702_clean.fastq" \
    -o "${ltde}/data/unicycler_assemblies/KBS0702_2"


unicycler \
    -1 "${ltde}/illumina_data/KBS0705_2015_SoilGenomes/KBS0705_R1_clean.fastq" \
    -2 "${ltde}/illumina_data/KBS0705_2015_SoilGenomes/KBS0705_R2_clean.fastq" \
    -l "${ltde}/data/nanopore_basecalled_bc_merged/KBS0705_clean.fastq" \
    -o "${ltde}/data/unicycler_assemblies/KBS0705"


unicycler \
    -1 "${ltde}/illumina_data/KBS0706_2015_SoilGenomes/KBS0706_R1_clean.fastq" \
    -2 "${ltde}/illumina_data/KBS0706_2015_SoilGenomes/KBS0706_R2_clean.fastq" \
    -l "${ltde}/data/nanopore_basecalled_bc_merged/KBS0706_clean.fastq" \
    -o "${ltde}/data/unicycler_assemblies/KBS0706"
