#!/bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=8,vmem=50gb,walltime=1:00:00
#PBS -M wrshoema@indiana.edu
#PBS -m abe
#PBS -j oe

module load bwa
module load samtools
module load bedtools
module load spades
module load python

# get RNA read
KBS0710="/N/dc2/projects/muri2/Task2/LTDE/data/align/KBS0710_NR_024911.fa"
KBS0721="/N/dc2/projects/muri2/Task2/LTDE/data/align/KBS0721_NR_114994.fa"

KBS0710_fastq="/N/dc2/projects/muri2/Task2/LTDE/2015_SoilGenomes/KBS0710/KBS0710.trim.interleaved.trim.fastq.pe.keep.pe"
KBS0721_fastq="/N/dc2/projects/muri2/Task2/LTDE/2015_SoilGenomes/KBS0721/KBS0721.trim.interleaved.trim.fastq.pe.keep.pe"

KBS0710_sam="/N/dc2/projects/muri2/Task2/LTDE/data/align/KBS0710_NR_024911.sam"
KBS0721_sam="/N/dc2/projects/muri2/Task2/LTDE/data/align/KBS0721_NR_114994.sam"

KBS0710_bam="/N/dc2/projects/muri2/Task2/LTDE/data/align/KBS0710_NR_024911.bam"
KBS0721_bam="/N/dc2/projects/muri2/Task2/LTDE/data/align/KBS0721_NR_114994.bam"

bwa index $KBS0710
bwa index $KBS0721
samtools faidx $KBS0710
samtools faidx $KBS0721



#KBS0710_bam_R1_fastq="/N/dc2/projects/muri2/Task2/LTDE/data/align/KBS0710_NR_024911_R1.fastq"
#KBS0710_bam_R2_fastq="/N/dc2/projects/muri2/Task2/LTDE/data/align/KBS0710_NR_024911_R2.fastq"
#KBS0721_bam_R1_fastq="/N/dc2/projects/muri2/Task2/LTDE/data/align/KBS0721_NR_114994_R1.fastq"
#KBS0721_bam_R2_fastq="/N/dc2/projects/muri2/Task2/LTDE/data/align/KBS0721_NR_114994_R2.fastq"

# for a single interleaved reads
# -p = paired end alignment
# -q phred quality score
# -F reads that mapped to the reference

samtools index $KBS0710_bam
samtools index $KBS0721_bam

bwa mem -M -t 4 -p $KBS0710 $KBS0710_fastq | samtools view -q 30 -F 4 -bT $KBS0710 - \
    | samtools sort -o $KBS0710_bam
bwa mem -M -t 4 -p $KBS0721 $KBS0721_fastq | samtools view -q 30 -F 4 -bT $KBS0721 - \
    | samtools sort -o $KBS0721_bam

KBS0710_pileup="/N/dc2/projects/muri2/Task2/LTDE/data/align/KBS0710_NR_024911.pileup"
KBS0721_pileup="/N/dc2/projects/muri2/Task2/LTDE/data/align/KBS0721_NR_114994.pileup"

samtools mpileup -q 30 -f $KBS0710 -o $KBS0710_pileup $KBS0710_bam
samtools mpileup -q 30 -f $KBS0721 -o $KBS0721_pileup $KBS0721_bam


#bedtools bamtofastq -i $KBS0710_bam \
#                      -fq $KBS0710_bam_R1_fastq \
#                      -fq2 $KBS0710_bam_R2_fastq

#bedtools bamtofastq -i $KBS0721_bam \
#                      -fq $KBS0721_bam_R1_fastq \
#                      -fq2 $KBS0721_bam_R2_fastq


# assemble alignment....
#KBS0710_align="/N/dc2/projects/muri2/Task2/LTDE/data/align/KBS0710_NR_024911_spades"
#KBS0721_align="/N/dc2/projects/muri2/Task2/LTDE/data/align/KBS0721_NR_114994_spades"
#spades.py --careful -1 $KBS0710_bam_R1_fastq -2 $KBS0710_bam_R2_fastq \
#    -o $KBS0710_align
#spades.py --careful -1 $KBS0721_bam_R1_fastq -2 $KBS0721_bam_R2_fastq \
#    -o $KBS0721_align

#samtools tview $KBS0710_bam $KBS0710

#samtools tview $KBS0721_bam $KBS0721
