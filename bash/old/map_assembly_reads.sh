#!/bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=8,vmem=100gb,walltime=6:00:00
#PBS -M wrshoema@iu.edu
#PBS -m abe
#PBS -j oe


module load samtools
module load bwa

# concatenate human genome, phg:phiX174, rattus norvegicus, Rainbow trout, and yeast
# GCF_002163495.1_Omyk_1.0_genomic.fna
# GCF_000146045.2_R64_genomic.fna
# hg38.fna
# NC_001422.fna
# GCF_000001895.5_Rnor_6.0_genomic.fna

MINIMAP2=/N/u/wrshoema/Carbonate/minimap2-2.17_x64-linux/minimap2
#REF=/N/dc2/projects/muri2/Task2/LTDE/data/contaim_genomes/contam.fasta
#cat /N/dc2/projects/muri2/Task2/LTDE/data/contaim_genomes/*.fna > $REF

RAT=/N/dc2/projects/muri2/Task2/LTDE/data/contaim_genomes/GCF_000001895.5_Rnor_6.0_genomic.fna
TROUT=/N/dc2/projects/muri2/Task2/LTDE/data/contaim_genomes/GCF_002163495.1_Omyk_1.0_genomic.fna

samtools faidx $RAT
samtools faidx $TROUT
bwa index $RAT
bwa index $TROUT

nano_path=/N/dc2/projects/muri2/Task2/LTDE/data/nanopore_basecalled_bc_merged
ill_path=/N/dc2/projects/muri2/Task2/LTDE/illumina_data

$MINIMAP2 -ax map-ont $RAT "${nano_path}/KBS0706_clean.fastq" | samtools fastq -n -f 4 - > "${nano_path}/KBS0706_clean_noContam.fastq"
$MINIMAP2 -ax map-ont $RAT "${nano_path}/KBS0713_clean.fastq" | samtools fastq -n -f 4 - > "${nano_path}/KBS0713_clean_noContam.fastq"
$MINIMAP2 -ax map-ont $RAT "${nano_path}/KBS0715_clean.fastq" | samtools fastq -n -f 4 - > "${nano_path}/KBS0715_clean_noContam.fastq"
$MINIMAP2 -ax map-ont $RAT "${nano_path}/KBS0802_clean.fastq" | samtools fastq -n -f 4 - > "${nano_path}/KBS0802_clean_noContam.fastq"
$MINIMAP2 -ax map-ont $TROUT "${nano_path}/KBS0722_clean.fastq" | samtools fastq -n -f 4 - > "${nano_path}/KBS0722_clean_noContam.fastq"

bwa mem $RAT "${ill_path}/KBS0706_2015_SoilGenomes/KBS0706_R1_clean.fastq" "${ill_path}/KBS0706_2015_SoilGenomes/KBS0706_R2_clean.fastq" \
    | samtools fastq -n -f 4 \
    -1 "${ill_path}/KBS0706_2015_SoilGenomes/KBS0706_R1_clean_noContam.fastq" \
    -2 "${ill_path}/KBS0706_2015_SoilGenomes/KBS0706_R2_clean_noContam.fastq"


bwa mem $RAT "${ill_path}/KBS0713_2015_SoilGenomes/KBS0713_R1_clean.fastq" "${ill_path}/KBS0713_2015_SoilGenomes/KBS0713_R2_clean.fastq" \
    | samtools fastq -n -f 4 \
    -1 "${ill_path}/KBS0713_2015_SoilGenomes/KBS0713_R1_clean_noContam.fastq" \
    -2 "${ill_path}/KBS0713_2015_SoilGenomes/KBS0713_R2_clean_noContam.fastq"


bwa mem $RAT "${ill_path}/KBS0715_2015_SoilGenomes/KBS0715_R1_clean.fastq" "${ill_path}/KBS0715_2015_SoilGenomes/KBS0715_R2_clean.fastq" \
    | samtools fastq -n -f 4 \
    -1 "${ill_path}/KBS0715_2015_SoilGenomes/KBS0715_R1_clean_noContam.fastq" \
    -2 "${ill_path}/KBS0715_2015_SoilGenomes/KBS0715_R2_clean_noContam.fastq"


bwa mem $RAT "${ill_path}/KBS0802_2015_SoilGenomes/KBS0802_R1_clean.fastq" "${ill_path}/KBS0802_2015_SoilGenomes/KBS0802_R2_clean.fastq" \
    | samtools fastq -n -f 4 \
    -1 "${ill_path}/KBS0802_2015_SoilGenomes/KBS0802_R1_clean_noContam.fastq" \
    -2 "${ill_path}/KBS0802_2015_SoilGenomes/KBS0802_R2_clean_noContam.fastq"


bwa mem $TROUT "${ill_path}/KBS0722_2015_SoilGenomes/KBS0722_R1_clean.fastq" "${ill_path}/KBS0722_2015_SoilGenomes/KBS0722_R2_clean.fastq" \
    | samtools fastq -n -f 4 \
    -1 "${ill_path}/KBS0722_2015_SoilGenomes/KBS0722_R1_clean_noContam.fastq" \
    -2 "${ill_path}/KBS0722_2015_SoilGenomes/KBS0722_R2_clean_noContam.fastq"
