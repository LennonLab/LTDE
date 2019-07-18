#!/bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=8,vmem=100gb,walltime=4:00:00
#PBS -M wrshoema@iu.edu
#PBS -m abe
#PBS -j oe

module load unicycler

ltde=/N/dc2/projects/muri2/Task2/LTDE

unicycler \
    -1 "${ltde}/illumina_data/GSF966-1-Arthro-6k_S1_R1_001_clean.fastq" \
    -1 "${ltde}/illumina_data/GSF966-1-Arthro-6k_S1_R2_001_clean.fastq" \
    -l "{ltde}/data/nanopore_basecalled_bc_merged/KBS0703_clean.fastq" \
    -o "${ltde}/data/unicycler_assemblies/KBS0703"


unicycler \
    -1 "${ltde}/illumina_data/ATCC13985_D400_100/Ev13483_TAAGGCGAATCT-AAGGAGTA_L002_R1_001_clean.fastq" \
    -2 "${ltde}/illumina_data/ATCC13985_D400_100/Ev13483_TAAGGCGAATCT-AAGGAGTA_L002_R2_001_clean.fastq" \
    -l "${ltde}/data/nanopore_basecalled_bc_merged/ATCC13985_clean.fastq" \
    -o "${ltde}/data/unicycler_assemblies/ATCC13985"


unicycler \
    -1 "${ltde}/illumina_data/ATCC43928_D400_100/Ev43828_TAAGGCGAATCT-ACTGCATA_L002_R1_001_clean.fastq.gz" \
    -2 "${ltde}/illumina_data/ATCC43928_D400_100/Ev43828_TAAGGCGAATCT-ACTGCATA_L002_R2_001_clean.fastq.gz" \
    -l "${ltde}/data/nanopore_basecalled_bc_merged/ATCC43928_clean.fastq" \
    -o "${ltde}/data/unicycler_assemblies/ATCC43928"


unicycler \
    -1 "${ltde}/illumina_data/KBS0701_2015_SoilGenomes/KBS0701_R1_clean.fastq" \
    -2 "${ltde}/illumina_data/KBS0701_2015_SoilGenomes/KBS0701_R2_clean.fastq" \
    -l "${ltde}/data/nanopore_basecalled_bc_merged/KBS0701_clean.fastq" \
    -o "${ltde}/data/unicycler_assemblies/KBS0701"


unicycler \
    -1 "${ltde}/illumina_data/KBS0702_D400_100/Ev702_TAAGGCGAATCT-TAGATCGC_L002_R1_001_clean.fastq.gz" \
    -2 "${ltde}/illumina_data/KBS0702_D400_100/Ev702_TAAGGCGAATCT-TAGATCGC_L002_R2_001_clean.fastq.gz" \
    -l "${ltde}/data/nanopore_basecalled_bc_merged/KBS0702_clean.fastq" \
    -o "${ltde}/data/unicycler_assemblies/KBS0702"


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
    -o "${ltde}/data/unicycler_assemblies/KBS0711"


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


unicycler \
    -1 "${ltde}/illumina_data/KBS0721_2015_SoilGenomes/KBS0721_R1_clean.fastq" \
    -2 "${ltde}/illumina_data/KBS0721_2015_SoilGenomes/KBS0721_R2_clean.fastq" \
    -l "${ltde}/data/nanopore_basecalled_bc_merged/KBS0721_clean.fastq" \
    -o "${ltde}/data/unicycler_assemblies/KBS0721"


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
    -o "${ltde}/data/unicycler_assemblies/KBS0801"


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
