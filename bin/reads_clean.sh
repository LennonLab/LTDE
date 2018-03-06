#!/bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=8,vmem=100gb,walltime=10:00:00
#PBS -M wrshoema@indiana.edu
#PBS -m abe
#PBS -j oe

module load fastqc
module load cutadapt
module load python
module load java

# fastqc of raw data
mkdir -p /N/dc2/projects/muri2/Task2/LTDE/data/reads_raw_quality/

for file in "/N/dc2/projects/muri2/Task2/LTDE/data/reads_raw/"*.fastq.gz
do
  sample="$(echo "$file" | cut -d "/" -f10-10 | cut -d "_" -f1-1 | cut -d "-" -f2-3)"
  OUTdir="/N/dc2/projects/muri2/Task2/LTDE/data/reads_raw_quality/${sample}"
  mkdir -p $OUTdir
  fastqc "$file" --outdir=$OUTdir
done


# rename files to include adaptors
mkdir -p /N/dc2/projects/muri2/Task2/LTDE/data/reads_raw_rename/
#python /N/dc2/projects/muri2/Task2/LTDE/bin/rename_libraries.py

# trim data and remove adaptors.
mkdir -p "/N/dc2/projects/muri2/Task2/LTDE/data/reads_clean"

#for R1 in "/N/dc2/projects/muri2/Task2/LTDE/data/reads_raw_rename/"*"_R1_"*.fastq.gz
#do
#  line="$(  echo "$R1" | cut -d"/" -f10-10 | cut -d "_" -f1-1)"
#  mkdir -p "/N/dc2/projects/muri2/Task2/LTDE/data/reads_clean/${line}"
#  R2="/N/dc2/projects/muri2/Task2/LTDE/data/reads_raw_rename/${line}"*"_R2_"*
#  adaptor1="$( echo "${R1}" | cut -d"_" -f7-7 | cut -d"." -f1-1 | cut -d"-" -f1-1 )"
#  adaptor2="$( echo "${R1}" | cut -d"_" -f7-7 | cut -d"." -f1-1 | cut -d"-" -f2-2 )"
#  R1_name="$( echo "${R1}" | cut -d "." -f1-1 )"
#  R2_name="$( echo "${R2}" | cut -d "." -f1-1 )"
#  OutR1Paired="/N/dc2/projects/muri2/Task2/LTDE/data/reads_clean/${line}/${line}_clean_R1_paired.fastq.gz"
#  OutR2Paired="/N/dc2/projects/muri2/Task2/LTDE/data/reads_clean/${line}/${line}_clean_R2_paired.fastq.gz"
#  OutR1UnPaired="/N/dc2/projects/muri2/Task2/LTDE/data/reads_clean/${line}/${line}_clean_R1_unpaired.fastq.gz"
#  OutR2UnPaired="/N/dc2/projects/muri2/Task2/LTDE/data/reads_clean/${line}/${line}_clean_R2_unpaired.fastq.gz"
#  java -jar /N/dc2/projects/MicroEukMA/softwares/Trimmomatic-0.32/trimmomatic-0.32.jar \
#    PE -threads 4 $R1 $R2 $OutR1Paired $OutR1UnPaired $OutR2Paired $OutR2UnPaired \
#    ILLUMINACLIP:/N/dc2/projects/muri2/Task2/PoolPopSeq/bin/transposase.fa:2:30:10 \
#    LEADING:4 TRAILING:4 MINLEN:40 HEADCROP:15
#done


mkdir -p /N/dc2/projects/muri2/Task2/LTDE/data/reads_clean_quality/

for file in "/N/dc2/projects/muri2/Task2/LTDE/data/reads_clean/"*"/"*_paired.fastq.gz
do
  sample="$(echo "$file" | cut -d "/" -f10-10 )"
  OUTdir="/N/dc2/projects/muri2/Task2/LTDE/data/reads_clean_quality/${sample}"
  mkdir -p $OUTdir
  fastqc "$file" --outdir=$OUTdir
done
