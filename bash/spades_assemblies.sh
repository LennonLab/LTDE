#!/bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=8,vmem=100gb,walltime=60:00:00
#PBS -M wrshoema@iu.edu
#PBS -m abe
#PBS -j oe

module load spades

ltde=/N/dc2/projects/muri2/Task2/LTDE

spades.py --careful --cov-cutoff auto \
      --nanopore "${ltde}/data/nanopore_basecalled_bc_merged/ATCC13985_clean.fastq" \
      --pe1-1 "${ltde}/illumina_data/ATCC13985_D400_100/Ev13483_TAAGGCGAATCT-AAGGAGTA_L002_R1_001_clean.fastq.gz" \
      --pe1-2 "${ltde}/illumina_data/ATCC13985_D400_100/Ev13483_TAAGGCGAATCT-AAGGAGTA_L002_R2_001_clean.fastq.gz" \
      -o "${ltde}/data/spades_assemblies/ATCC13985"

spades.py --careful --cov-cutoff auto \
      --nanopore "${ltde}/data/nanopore_basecalled_bc_merged/ATCC43928_clean.fastq" \
      --pe1-1 "${ltde}/illumina_data/ATCC43928_D400_100/Ev43828_TAAGGCGAATCT-ACTGCATA_L002_R1_001_clean.fastq.gz" \
      --pe1-2 "${ltde}/illumina_data/ATCC43928_D400_100/Ev43828_TAAGGCGAATCT-ACTGCATA_L002_R2_001_clean.fastq.gz" \
      -o "${ltde}/data/spades_assemblies/ATCC43928"

spades.py --careful --cov-cutoff auto \
      --nanopore "${ltde}/data/nanopore_basecalled_bc_merged/KBS0701_clean.fastq" \
      --pe1-1 "${ltde}/illumina_data/KBS0701_2015_SoilGenomes/KBS0701_R1_clean.fastq" \
      --pe1-2 "${ltde}/illumina_data/KBS0701_2015_SoilGenomes/KBS0701_R2_clean.fastq" \
      -o "${ltde}/data/spades_assemblies/KBS0701"

spades.py --careful --cov-cutoff auto \
      --nanopore "${ltde}/data/nanopore_basecalled_bc_merged/KBS0702_clean.fastq" \
      --pe1-1 "${ltde}/illumina_data/KBS0702_2015_SoilGenomes/KBS0702_R1_clean.fastq" \
      --pe1-2 "${ltde}/illumina_data/KBS0702_2015_SoilGenomes/KBS0702_R2_clean.fastq" \
      --pe2-1 "${ltde}/illumina_data/KBS0702_D400_100/Ev702_TAAGGCGAATCT-TAGATCGC_L002_R1_001_clean.fastq.gz" \
      --pe2-2 "${ltde}/illumina_data/KBS0702_D400_100/Ev702_TAAGGCGAATCT-TAGATCGC_L002_R2_001_clean.fastq.gz" \
      -o "${ltde}/data/spades_assemblies/KBS0702"

#spades.py --careful --cov-cutoff auto \
#      --nanopore "${ltde}/data/nanopore_basecalled_bc_merged/KBS0703_clean.fastq" \
#      --pe1-1 "${ltde}/illumina_data/KBS0703_GSF911/GSF911-Ar_S10_L001_R1_001_clean.fastq" \
#      --pe1-2 "${ltde}/illumina_data/KBS0703_GSF911/GSF911-Ar_S10_L001_R1_001_clean.fastq" \
#      -o "${ltde}/data/spades_assemblies/KBS0703"

spades.py --careful --cov-cutoff auto \
      --nanopore "${ltde}/data/nanopore_basecalled_bc_merged/KBS0705_clean.fastq" \
      --pe1-1 "${ltde}/illumina_data/KBS0705_2015_SoilGenomes/KBS0705_R1_clean.fastq" \
      --pe1-2 "${ltde}/illumina_data/KBS0705_2015_SoilGenomes/KBS0705_R2_clean.fastq" \
      -o "${ltde}/data/spades_assemblies/KBS0705"

spades.py --careful --cov-cutoff auto \
      --nanopore "${ltde}/data/nanopore_basecalled_bc_merged/KBS0706_clean.fastq" \
      --pe1-1 "${ltde}/illumina_data/KBS0706_2015_SoilGenomes/KBS0706_R1_clean.fastq" \
      --pe1-2 "${ltde}/illumina_data/KBS0706_2015_SoilGenomes/KBS0706_R2_clean.fastq" \
      -o "${ltde}/data/spades_assemblies/KBS0706"

spades.py --careful --cov-cutoff auto \
      --nanopore "${ltde}/data/nanopore_basecalled_bc_merged/KBS0707_clean.fastq" \
      --pe1-1 "${ltde}/illumina_data/KBS0707_D400_100/Ev707_TAAGGCGAATCT-CTCTCTAT_L002_R1_001_clean.fastq.gz" \
      --pe1-2 "${ltde}/illumina_data/KBS0707_D400_100/Ev707_TAAGGCGAATCT-CTCTCTAT_L002_R2_001_clean.fastq.gz" \
      -o "${ltde}/data/spades_assemblies/KBS0707"

spades.py --careful --cov-cutoff auto \
      --nanopore "${ltde}/data/nanopore_basecalled_bc_merged/KBS0710_clean.fastq" \
      --pe1-1 "${ltde}/illumina_data/KBS0710_2015_SoilGenomes/KBS0710_R1_clean.fastq" \
      --pe1-2 "${ltde}/illumina_data/KBS0710_2015_SoilGenomes/KBS0710_R2_clean.fastq" \
      -o "${ltde}/data/spades_assemblies/KBS0710"

spades.py --careful --cov-cutoff auto \
      --nanopore "${ltde}/data/nanopore_basecalled_bc_merged/KBS0711_clean.fastq" \
      --pe1-1 "${ltde}/illumina_data/KBS0711_2015_SoilGenomes/KBS0711_R1_clean.fastq" \
      --pe1-2 "${ltde}/illumina_data/KBS0711_2015_SoilGenomes/KBS0711_R2_clean.fastq" \
      --pe2-1 "${ltde}/illumina_data/KBS0711_GSF911/GSF911-711_S1_L001_R1_001_clean.fastq" \
      --pe2-2 "${ltde}/illumina_data/KBS0711_GSF911/GSF911-711_S1_L001_R2_001_clean.fastq" \
      -o "${ltde}/data/spades_assemblies/KBS0711"

spades.py --careful --cov-cutoff auto \
      --nanopore "${ltde}/data/nanopore_basecalled_bc_merged/KBS0712_clean.fastq" \
      --pe1-1 "${ltde}/illumina_data/KBS0712_D400_100/Ev712_TAAGGCGAATCT-TATCCTCT_L002_R1_001_clean.fastq.gz" \
      --pe1-2 "${ltde}/illumina_data/KBS0712_D400_100/Ev712_TAAGGCGAATCT-TATCCTCT_L002_R2_001_clean.fastq.gz" \
      -o "${ltde}/data/spades_assemblies/KBS0712"

spades.py --careful --cov-cutoff auto \
      --nanopore "${ltde}/data/nanopore_basecalled_bc_merged/KBS0713_clean.fastq" \
      --pe1-1 "${ltde}/illumina_data/KBS0713_2015_SoilGenomes/KBS0713_R1_clean.fastq" \
      --pe1-2 "${ltde}/illumina_data/KBS0713_2015_SoilGenomes/KBS0713_R2_clean.fastq" \
      -o "${ltde}/data/spades_assemblies/KBS0713"

spades.py --careful --cov-cutoff auto \
      --nanopore "${ltde}/data/nanopore_basecalled_bc_merged/KBS0714_clean.fastq" \
      --pe1-1 "${ltde}/illumina_data/KBS0714_2015_SoilGenomes/KBS0714_R1_clean.fastq" \
      --pe1-2 "${ltde}/illumina_data/KBS0714_2015_SoilGenomes/KBS0714_R2_clean.fastq" \
      -o "${ltde}/data/spades_assemblies/KBS0714"


spades.py --careful --cov-cutoff auto \
      --nanopore "${ltde}/data/nanopore_basecalled_bc_merged/KBS0715_clean.fastq" \
      --pe1-1 "${ltde}/illumina_data/KBS0715_2015_SoilGenomes/KBS0715_R1_clean.fastq" \
      --pe1-2 "${ltde}/illumina_data/KBS0715_2015_SoilGenomes/KBS0715_R2_clean.fastq" \
      -o "${ltde}/data/spades_assemblies/KBS0715"

spades.py --careful --cov-cutoff auto \
      --nanopore "${ltde}/data/nanopore_basecalled_bc_merged/KBS0721_clean.fastq" \
      --pe1-1 "${ltde}/illumina_data/KBS0721_2015_SoilGenomes/KBS0721_R1_clean.fastq" \
      --pe1-2 "${ltde}/illumina_data/KBS0721_2015_SoilGenomes/KBS0721_R2_clean.fastq" \
      --pe2-1 "${ltde}/illumina_data/KBS0721_GSF911/GSF911-Fl_S9_L001_R1_001_clean.fastq" \
      --pe2-2 "${ltde}/illumina_data/KBS0721_GSF911/GSF911-Fl_S9_L001_R2_001_clean.fastq" \
      -o "${ltde}/data/spades_assemblies/KBS0721"

spades.py --careful --cov-cutoff auto \
      --nanopore "${ltde}/data/nanopore_basecalled_bc_merged/KBS0722_clean.fastq" \
      --pe1-1 "${ltde}/illumina_data/KBS0722_2015_SoilGenomes/KBS0722_R1_clean.fastq" \
      --pe1-2 "${ltde}/illumina_data/KBS0722_2015_SoilGenomes/KBS0722_R2_clean.fastq" \
      -o "${ltde}/data/spades_assemblies/KBS0722"

spades.py --careful --cov-cutoff auto \
      --nanopore "${ltde}/data/nanopore_basecalled_bc_merged/KBS0724_clean.fastq" \
      --pe1-1 "${ltde}/illumina_data/KBS0724_2015_SoilGenomes/KBS0724_R1_clean.fastq" \
      --pe1-2 "${ltde}/illumina_data/KBS0724_2015_SoilGenomes/KBS0724_R2_clean.fastq" \
      -o "${ltde}/data/spades_assemblies/KBS0724"

spades.py --careful --cov-cutoff auto \
      --nanopore "${ltde}/data/nanopore_basecalled_bc_merged/KBS0725_clean.fastq" \
      --pe1-1 "${ltde}/illumina_data/KBS0725_2015_SoilGenomes/KBS0725_R1_clean.fastq" \
      --pe1-2 "${ltde}/illumina_data/KBS0725_2015_SoilGenomes/KBS0725_R2_clean.fastq" \
      -o "${ltde}/data/spades_assemblies/KBS0725"

spades.py --careful --cov-cutoff auto \
      --nanopore "${ltde}/data/nanopore_basecalled_bc_merged/KBS0727_clean.fastq" \
      --pe1-1 "${ltde}/illumina_data/KBS0727_2015_SoilGenomes/KBS0727_R1_clean.fastq" \
      --pe1-2 "${ltde}/illumina_data/KBS0727_2015_SoilGenomes/KBS0727_R2_clean.fastq" \
      -o "${ltde}/data/spades_assemblies/KBS0727"

spades.py --careful --cov-cutoff auto \
      --nanopore "${ltde}/data/nanopore_basecalled_bc_merged/KBS0801_clean.fastq" \
      --pe1-1 "${ltde}/illumina_data/KBS0801_B_2015_SoilGenomes/KBS0801B_R1_clean.fastq" \
      --pe1-2 "${ltde}/illumina_data/KBS0801_B_2015_SoilGenomes/KBS0801B_R2_clean.fastq" \
      --pe2-1 "${ltde}/illumina_data/KBS0801_D400_100/Ev801_TAAGGCGAATCT-AGAGTAGA_L002_R1_001_clean.fastq.gz" \
      --pe2-2 "${ltde}/illumina_data/KBS0801_D400_100/Ev801_TAAGGCGAATCT-AGAGTAGA_L002_R2_001_clean.fastq.gz" \
      -o "${ltde}/data/spades_assemblies/KBS0801"

spades.py --careful --cov-cutoff auto \
      --nanopore "${ltde}/data/nanopore_basecalled_bc_merged/KBS0802_clean.fastq" \
      --pe1-1 "${ltde}/illumina_data/KBS0802_2015_SoilGenomes/KBS0802_R1_clean.fastq" \
      --pe1-2 "${ltde}/illumina_data/KBS0802_2015_SoilGenomes/KBS0802_R2_clean.fastq" \
      -o "${ltde}/data/spades_assemblies/KBS0802"

spades.py --careful --cov-cutoff auto \
      --nanopore "${ltde}/data/nanopore_basecalled_bc_merged/KBS0812_clean.fastq" \
      --pe1-1 "${ltde}/illumina_data/KBS0812_D400_100/Ev812_TAAGGCGAATCT-GTAAGGAG_L002_R1_001_clean.fastq.gz" \
      --pe1-2 "${ltde}/illumina_data/KBS0812_D400_100/Ev812_TAAGGCGAATCT-GTAAGGAG_L002_R2_001_clean.fastq.gz" \
      -o "${ltde}/data/spades_assemblies/KBS0812"


# just get contigs
for d in /N/dc2/projects/muri2/Task2/LTDE/data/spades_assemblies/* ; do
    #echo "${d}/contigs.fasta"
    strain="$(echo "$d" | cut -d "/" -f10-10 )"
    cp "${d}/contigs.fasta" "/N/dc2/projects/muri2/Task2/LTDE/data/spades_assemblies_contigs/${strain}.fasta"
done


for d in /N/dc2/projects/muri2/Task2/LTDE/data/spades_assemblies/* ; do
    strain="$(echo "$d" | cut -d "/" -f10-10 )"
    cp "${d}/contigs.fasta" "/N/dc2/projects/muri2/Task2/LTDE/data/spades_assemblies_contigs/${strain}.fasta"
done


for d in /Users/WRShoemaker/GitHub/LTDE/data/genomes/nanopore_hybrid/*.fasta; do
    strain="$(echo "$d" | cut -d "/" -f9-9 | cut -d "." -f1-1 )"
    prokka --compliant --centre IUB --outdir PRJNA539822 --locustag FCE86 \
            --prefix FCE86-Genome --outdir "/Users/WRShoemaker/GitHub/LTDE/data/genomes/nanopore_hybrid_annotated/${strain}"\
            $d
done


for d in /Users/WRShoemaker/GitHub/LTDE/data/genomes/nanopore_hybrid_annotated/*; do
    strain="$(echo "$d" | cut -d "/" -f9-9 )"
    echo $strain
    REF_OUT="~/LTDE/data/genomes/nanopore_hybrid_annotated_cogs/${strain}_reformat"
    anvi-script-reformat-fasta $d/FCE86-Genome.fna -o $REF_OUT.fna -l 0 --simplify-names
    anvi-gen-contigs-database -f $REF_OUT.fna -o $REF_OUT.db
    anvi-run-ncbi-cogs -c $REF_OUT.db --num-threads 20
    anvi-export-functions -c $REF_OUT.db -o $REF_OUT.txt
done
