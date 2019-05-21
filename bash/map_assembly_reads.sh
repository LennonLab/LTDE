#!/bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=8,vmem=100gb,walltime=4:00:00
#PBS -M wrshoema@iu.edu
#PBS -m abe
#PBS -j oe


module load samtools

# concatenate human genome, phg:phiX174, rattus norvegicus, Rainbow trout, and yeast
# GCF_002163495.1_Omyk_1.0_genomic.fna
# GCF_000146045.2_R64_genomic.fna
# hg38.fna
# NC_001422.fna
# GCF_000001895.5_Rnor_6.0_genomic.fna

MINIMAP2=/N/u/wrshoema/Carbonate/minimap2-2.17_x64-linux/minimap2
REF=/N/dc2/projects/muri2/Task2/LTDE/data/contaim_genomes/contam.fasta
cat /N/dc2/projects/muri2/Task2/LTDE/data/contaim_genomes/*.fna > $REF

samtools faidx $REF

declare -a to_clean=("KBS0701" "KBS0701" "KBS0703" "KBS0705" "ATCC13985")

for i in "${to_clean[@]}"
do
    FASTQ="/N/dc2/projects/muri2/Task2/LTDE/data/nanopore_basecalled_pcbc_rn/${i}.fastq"
    FASTQ_CLEAN="/N/dc2/projects/muri2/Task2/LTDE/data/nanopore_basecalled_pcbc_rn/${i}_clean.fastq"
    $MINIMAP2 -ax map-ont $REF $FASTQ | samtools fastq -n -f 4 - > $FASTQ_CLEAN
done
