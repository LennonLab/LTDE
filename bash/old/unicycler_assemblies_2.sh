#!/bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=8,vmem=100gb,walltime=120:00:00
#PBS -M wrshoema@iu.edu
#PBS -m abe
#PBS -j oe

module load unicycler

ltde=/N/dc2/projects/muri2/Task2/LTDE



unicycler \
    -1 "${ltde}/illumina_data/KBS0707_D400_100/Ev707_TAAGGCGAATCT-CTCTCTAT_L002_R1_001_clean.fastq.gz" \
    -2 "${ltde}/illumina_data/KBS0707_D400_100/Ev707_TAAGGCGAATCT-CTCTCTAT_L002_R2_001_clean.fastq.gz" \
    -l "${ltde}/data/nanopore_basecalled_bc_merged/KBS0707_clean.fastq" \
    -o "${ltde}/data/unicycler_assemblies/KBS0707"


unicycler \
    -1 "${ltde}/illumina_data/KBS0710_2015_SoilGenomes/KBS0710_R1_clean.fastq" \
    -2 "${ltde}/illumina_data/KBS0710_2015_SoilGenomes/KBS0710_R2_clean.fastq" \
    -l "${ltde}/data/nanopore_basecalled_bc_merged/KBS0710_clean.fastq" \
    -o "${ltde}/data/unicycler_assemblies/KBS0710"


unicycler \
    -1 "${ltde}/illumina_data/KBS0711_2015_SoilGenomes/KBS0711_R1_clean.fastq" \
    -2 "${ltde}/illumina_data/KBS0711_2015_SoilGenomes/KBS0711_R2_clean.fastq" \
    -l "${ltde}/data/nanopore_basecalled_bc_merged/KBS0711_clean.fastq" \
    -o "${ltde}/data/unicycler_assemblies/KBS0711_1"


unicycler \
    -1 "${ltde}/illumina_data/KBS0711_GSF911/GSF911-711_S1_L001_R1_001_clean.fastq" \
    -2 "${ltde}/illumina_data/KBS0711_GSF911/GSF911-711_S1_L001_R2_001_clean.fastq" \
    -l "${ltde}/data/nanopore_basecalled_bc_merged/KBS0711_clean.fastq" \
    -o "${ltde}/data/unicycler_assemblies/KBS0711_2"




unicycler \
    -1 "${ltde}/illumina_data/KBS0712_D400_100/Ev712_TAAGGCGAATCT-TATCCTCT_L002_R1_001_clean.fastq" \
    -2 "${ltde}/illumina_data/KBS0712_D400_100/Ev712_TAAGGCGAATCT-TATCCTCT_L002_R2_001_clean.fastq" \
    -l "${ltde}/data/nanopore_basecalled_bc_merged/KBS0712_clean.fastq" \
    -o "${ltde}/data/unicycler_assemblies/KBS0712"


unicycler \
    -1 "${ltde}/illumina_data/KBS0713_2015_SoilGenomes/KBS0713_R1_clean.fastq" \
    -2 "${ltde}/illumina_data/KBS0713_2015_SoilGenomes/KBS0713_R2_clean.fastq" \
    -l "${ltde}/data/nanopore_basecalled_bc_merged/KBS0713_clean.fastq" \
    -o "${ltde}/data/unicycler_assemblies/KBS0713"


unicycler \
    -1 "${ltde}/illumina_data/KBS0714_2015_SoilGenomes/KBS0714_R1_clean.fastq" \
    -2 "${ltde}/illumina_data/KBS0714_2015_SoilGenomes/KBS0714_R2_clean.fastq" \
    -l "${ltde}/data/nanopore_basecalled_bc_merged/KBS0714_clean.fastq" \
    -o "${ltde}/data/unicycler_assemblies/KBS0714"


unicycler \
    -1 "${ltde}/illumina_data/KBS0715_2015_SoilGenomes/KBS0715_R1_clean.fastq" \
    -2 "${ltde}/illumina_data/KBS0715_2015_SoilGenomes/KBS0715_R2_clean.fastq" \
    -l "${ltde}/data/nanopore_basecalled_bc_merged/KBS0715_clean.fastq" \
    -o "${ltde}/data/unicycler_assemblies/KBS0715"
