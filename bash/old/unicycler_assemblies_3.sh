#!/bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=8,vmem=100gb,walltime=120:00:00
#PBS -M wrshoema@iu.edu
#PBS -m abe
#PBS -j oe

module load unicycler

ltde=/N/dc2/projects/muri2/Task2/LTDE



unicycler \
    -1 "${ltde}/illumina_data/KBS0721_2015_SoilGenomes/KBS0721_R1_clean.fastq" \
    -2 "${ltde}/illumina_data/KBS0721_2015_SoilGenomes/KBS0721_R2_clean.fastq" \
    -l "${ltde}/data/nanopore_basecalled_bc_merged/KBS0721_clean.fastq" \
    -o "${ltde}/data/unicycler_assemblies/KBS0721_1"


unicycler \
    -1 "${ltde}/illumina_data/KBS0721_GSF911/GSF911-Fl_S9_L001_R1_001_clean.fastq" \
    -2 "${ltde}/illumina_data/KBS0721_GSF911/GSF911-Fl_S9_L001_R2_001_clean.fastq" \
    -l "${ltde}/data/nanopore_basecalled_bc_merged/KBS0721_clean.fastq" \
    -o "${ltde}/data/unicycler_assemblies/KBS0721_2"


unicycler \
    -1 "${ltde}/illumina_data/KBS0722_2015_SoilGenomes/KBS0722_R1_clean.fastq" \
    -2 "${ltde}/illumina_data/KBS0722_2015_SoilGenomes/KBS0722_R2_clean.fastq" \
    -l "${ltde}/data/nanopore_basecalled_bc_merged/KBS0722_clean.fastq" \
    -o "${ltde}/data/unicycler_assemblies/KBS0722"


unicycler \
    -1 "${ltde}/illumina_data/KBS0724_2015_SoilGenomes/KBS0724_R1_clean.fastq" \
    -2 "${ltde}/illumina_data/KBS0724_2015_SoilGenomes/KBS0724_R2_clean.fastq" \
    -l "${ltde}/data/nanopore_basecalled_bc_merged/KBS0724_clean.fastq" \
    -o "${ltde}/data/unicycler_assemblies/KBS0724"


unicycler \
    -1 "${ltde}/illumina_data/KBS0725_2015_SoilGenomes/KBS0725_R1_clean.fastq" \
    -2 "${ltde}/illumina_data/KBS0725_2015_SoilGenomes/KBS0725_R2_clean.fastq" \
    -l "${ltde}/data/nanopore_basecalled_bc_merged/KBS0725_clean.fastq" \
    -o "${ltde}/data/unicycler_assemblies/KBS0725"


unicycler \
    -1 "${ltde}/illumina_data/KBS0727_2015_SoilGenomes/KBS0727_R1_clean.fastq" \
    -2 "${ltde}/illumina_data/KBS0727_2015_SoilGenomes/KBS0727_R2_clean.fastq" \
    -l "${ltde}/data/nanopore_basecalled_bc_merged/KBS0727_clean.fastq" \
    -o "${ltde}/data/unicycler_assemblies/KBS0727"


unicycler \
    -1 "${ltde}/illumina_data/KBS0801_B_2015_SoilGenomes/KBS0801B_R1_clean.fastq" \
    -2 "${ltde}/illumina_data/KBS0801_B_2015_SoilGenomes/KBS0801B_R2_clean.fastq" \
    -l "${ltde}/data/nanopore_basecalled_bc_merged/KBS0801_clean.fastq" \
    -o "${ltde}/data/unicycler_assemblies/KBS0801_1"


unicycler \
    -1 "${ltde}/illumina_data/KBS0801_D400_100/Ev801_TAAGGCGAATCT-AGAGTAGA_L002_R1_001_clean.fastq.gz" \
    -2 "${ltde}/illumina_data/KBS0801_D400_100/Ev801_TAAGGCGAATCT-AGAGTAGA_L002_R2_001_clean.fastq.gz" \
    -l "${ltde}/data/nanopore_basecalled_bc_merged/KBS0801_clean.fastq" \
    -o "${ltde}/data/unicycler_assemblies/KBS0801_2"


unicycler \
    -1 "${ltde}/illumina_data/KBS0802_2015_SoilGenomes/KBS0802_R1_clean.fastq" \
    -2 "${ltde}/illumina_data/KBS0802_2015_SoilGenomes/KBS0802_R2_clean.fastq" \
    -l "${ltde}/data/nanopore_basecalled_bc_merged/KBS0802_clean.fastq" \
    -o "${ltde}/data/unicycler_assemblies/KBS0802"


unicycler \
    -1 "${ltde}/illumina_data/KBS0812_D400_100/Ev812_TAAGGCGAATCT-GTAAGGAG_L002_R1_001_clean.fastq.gz" \
    -2 "${ltde}/illumina_data/KBS0812_D400_100/Ev812_TAAGGCGAATCT-GTAAGGAG_L002_R2_001_clean.fastq.gz" \
    -l "${ltde}/data/nanopore_basecalled_bc_merged/KBS0812_clean.fastq" \
    -o "${ltde}/data/unicycler_assemblies/KBS0812"
