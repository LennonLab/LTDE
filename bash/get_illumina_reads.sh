#!/bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=8,vmem=100gb,walltime=4:00:00
#PBS -M wrshoema@iu.edu
#PBS -m abe
#PBS -j oe

#module unload python
#module load python/3.6.1
module load python/2.7.13
module load fastqc
module load cutadapt/intel/1.16


for folder in /N/dc2/projects/muri2/Task2/LTDE/illumina_data/*/
do
  if [[ $folder == *"D400_100"* ]]; then
    sample_sheet="${folder}SampleSheet.csv"
    adapt1="$(sed '2q;d' ${sample_sheet} | cut -d"," -f5 | cut -d"-" -f1)"
    adapt2="$(sed '2q;d' ${sample_sheet} | cut -d"," -f5 | cut -d"-" -f2)"
    R1="${folder}"*"_R1_001.fastq.gz"
    R2="${folder}"*"_R2_001.fastq.gz"
    #zgrep ^@HSQ $R1 | cut -d":" -f10 | sort -r | uniq > "${folder}"barcodes_R1.txt
    #zgrep ^@HSQ $R2 | cut -d":" -f10 | sort -r | uniq > "${folder}"barcodes_R2.txt
    #cat "${folder}barcodes_R"?".txt" > "${folder}barcodes.txt"
    #sort -u "${folder}barcodes.txt" > "${folder}barcodes_unique.txt"
    #python /N/dc2/projects/muri2/Task2/LTDE/bash/add_fasta_headers.py "${folder}barcodes_unique.txt" "${folder}barcodes_unique.fa"

    R1_name="$( echo $R1 | cut -d "." -f1-1 )"
    R2_name="$( echo $R2 | cut -d "." -f1-1 )"

    if [[ $folder == *"ATCC13985_D400_100"* ]]; then
        #-g file:"${folder}barcodes_unique.fa" \
        cutadapt -q 30,30 -u 20 -a $adapt1 -A $adapt2 \
                  -a file:/N/dc2/projects/muri2/Task2/LTDE/bash/adapters_ATCC13985.fasta \
                  -A file:/N/dc2/projects/muri2/Task2/LTDE/bash/adapters_ATCC13985.fasta \
                  --minimum-length 20 \
                  -n 3 \
                  -o "${R1_name}_clean.fastq" \
                  -p "${R2_name}_clean.fastq" $R1 $R2

    elif [[ $folder == *"KBS0712_D400_100"* ]]; then
        cutadapt -q 30,30 -u 20 -a $adapt1 -A $adapt2 \
                  -a file:/N/dc2/projects/muri2/Task2/LTDE/bash/adapters_KBS0712.fasta \
                  -A file:/N/dc2/projects/muri2/Task2/LTDE/bash/adapters_KBS0712.fasta \
                  --minimum-length 20 \
                  -n 3 \
                  -o "${R1_name}_clean.fastq" \
                  -p "${R2_name}_clean.fastq" $R1 $R2

    else
        continue
        cutadapt -q 30,30 -u 20 -a $adapt1 -A $adapt2 \
                  --minimum-length 20 -o "${R1_name}_clean.fastq" \
                  -p "${R2_name}_clean.fastq" $R1 $R2
    fi


  elif [[ $folder == *"_2015_SoilGenomes"* ]]; then
    adapt1="AGATCGGAAGAGC"
    adapt2="AGATCGGAAGAGC"
    R1="${folder}"*"_R1.fastq"
    R2="${folder}"*"_R2.fastq"

    #grep ^@HWI $R1 | cut -d":" -f10 | sort -r | uniq > "${folder}"barcodes_R1.txt
    #grep ^@HWI $R2 | cut -d":" -f10 | sort -r | uniq > "${folder}"barcodes_R2.txt
    #cat "${folder}barcodes_R"?".txt" > "${folder}barcodes.txt"
    #sort -u "${folder}barcodes.txt" > "${folder}barcodes_unique.txt"
    #python /N/dc2/projects/muri2/Task2/LTDE/bash/add_fasta_headers.py "${folder}barcodes_unique.txt" "${folder}barcodes_unique.fa"

    R1_name="$( echo $R1 | cut -d "." -f1-1 )"
    R2_name="$( echo $R2 | cut -d "." -f1-1 )"

    #-g file:"${folder}barcodes_unique.fa" \
    if [[ $folder == *"KBS0721_2015_SoilGenomes"* ]]; then
        cutadapt -q 30,30 -u 20 -a $adapt1 -A $adapt2 \
                  -a file:/N/dc2/projects/muri2/Task2/LTDE/bash/adapters_KBS0721.fasta \
                  -A file:/N/dc2/projects/muri2/Task2/LTDE/bash/adapters_KBS0721.fasta \
                  --minimum-length 20 \
                  -n 3 \
                  -o "${R1_name}_clean.fastq" \
                  -p "${R2_name}_clean.fastq" $R1 $R2

    else
        continue
        cutadapt -q 30,30 -u 20 -a $adapt1 -A $adapt2 \
                  --minimum-length 20 -o "${R1_name}_clean.fastq" \
                  -p "${R2_name}_clean.fastq" $R1 $R2
    fi


  elif [[ $folder == *"_GSF911"* ]]; then
    R1="${folder}"*"_R1_001.fastq.gz"
    R2="${folder}"*"_R2_001.fastq.gz"

    # no barcode listed
    R1_name="$( echo $R1 | cut -d "." -f1-1 )"
    R2_name="$( echo $R2 | cut -d "." -f1-1 )"
    #adapt=/N/dc2/projects/muri2/Task2/LTDE/bash/GSF911_adapters.fasta

    if [[ $folder == *"KBS0703_GSF911"* ]]; then
        cutadapt -q 30,30 -u 20 \
                  -a file:/N/dc2/projects/muri2/Task2/LTDE/bash/GSF911_adapters.fasta \
                  -a file:/N/dc2/projects/muri2/Task2/LTDE/bash/adapters_KBS0703.fasta \
                  -A file:/N/dc2/projects/muri2/Task2/LTDE/bash/adapters_KBS0703.fasta \
                  --minimum-length 20 \
                  -n 3 \
                  -o "${R1_name}_clean.fastq" \
                  -p "${R2_name}_clean.fastq" $R1 $R2

    elif [[ $folder == *"KBS0721_GSF911"* ]]; then
        cutadapt -q 30,30 -u 20 \
                  -a file:/N/dc2/projects/muri2/Task2/LTDE/bash/GSF911_adapters.fasta \
                  -a file:/N/dc2/projects/muri2/Task2/LTDE/bash/adapters_KBS0721.fasta \
                  -A file:/N/dc2/projects/muri2/Task2/LTDE/bash/adapters_KBS0721.fasta \
                  --minimum-length 20 \
                  -n 3 \
                  -o "${R1_name}_clean.fastq" \
                  -p "${R2_name}_clean.fastq" $R1 $R2

    else
        continue
        cutadapt -q 30,30 -u 20 \
                  -a file:/N/dc2/projects/muri2/Task2/LTDE/bash/GSF911_adapters.fasta \
                  --minimum-length 20 -o "${R1_name}_clean.fastq" \
                  -p "${R2_name}_clean.fastq" $R1 $R2
    fi

  else
    continue

  #fastqc $R1 > "${folder}"R1_results.out
  #fastqc $R2 > "${folder}"R2_results.out

  #fastqc "${R1_name}_clean.fastq.gz" > "${folder}"R1_clean_results.out
  #fastqc "${R2_name}_clean.fastq.gz" > "${folder}"R2_clean_results.out
  fi
done
