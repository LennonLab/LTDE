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

declare -a strains=("ATCC13985" "ATCC43928" "KBS0701" "KBS0702" "KBS0703" \
                    "KBS0705" "KBS0706" "KBS0707" "KBS0710" "KBS0711" \
                    "KBS0712" "KBS0713" "KBS0714" "KBS0715" "KBS0721" \
                    "KBS0722" "KBS0724" "KBS0725" "KBS0727" "KBS0801" \
                    "KBS0802" "KBS0812")

declare -a libraries=("ATCC13985_D400_100" "ATCC43928_D400_100" \
                      "KBS0701_2015_SoilGenomes" "KBS0702_2015_SoilGenomes" \
                      "KBS0702_D400_100" "KBS0703_GSF911" "KBS0705_2015_SoilGenomes" \
                      "KBS0706_2015_SoilGenomes" "KBS0707_D400_100" "KBS0710_2015_SoilGenomes" \
                      "KBS0711_2015_SoilGenomes" "KBS0711_GSF911" "KBS0712_D400_100" \
                      "KBS0713_2015_SoilGenomes" "KBS0714_2015_SoilGenomes" \
                      "KBS0715_2015_SoilGenomes" "KBS0721_2015_SoilGenomes" \
                      "KBS0721_GSF911" "KBS0722_2015_SoilGenomes" "KBS0724_2015_SoilGenomes" \
                      "KBS0725_2015_SoilGenomes" "KBS0727_2015_SoilGenomes" \
                      "KBS0801_B_2015_SoilGenomes" "KBS0801_D400_100" \
                      "KBS0802_2015_SoilGenomes" "KBS0812_D400_100")

#for folder in /N/dc2/projects/muri2/Task2/LTDE/illumina_data/*/
#do
#  echo $folder
#done

for folder in /N/dc2/projects/muri2/Task2/LTDE/illumina_data/*/
do
  if [[ $folder == *"D400_100"* ]]; then
    continue
    sample_sheet="${folder}SampleSheet.csv"
    adapt1="$(sed '2q;d' ${sample_sheet} | cut -d"," -f5 | cut -d"-" -f1)"
    adapt2="$(sed '2q;d' ${sample_sheet} | cut -d"," -f5 | cut -d"-" -f2)"
    R1="${folder}"*"_R1_001.fastq.gz"
    R2="${folder}"*"_R2_001.fastq.gz"
    zgrep ^@HSQ $R1 | cut -d":" -f10 | sort -r | uniq > "${folder}"barcodes_R1.txt
    zgrep ^@HSQ $R2 | cut -d":" -f10 | sort -r | uniq > "${folder}"barcodes_R2.txt
    cat "${folder}barcodes_R"?".txt" > "${folder}barcodes.txt"
    #sort -u "${folder}barcodes.txt" >> "${folder}barcodes_unique.txt"
    sort "${folder}barcodes.txt" | uniq > "${folder}barcodes_unique.txt"
    python /N/dc2/projects/muri2/Task2/LTDE/bash/add_fasta_headers.py "${folder}barcodes_unique.txt" "${folder}barcodes_unique.fa"

    R1_name="$( echo $R1 | cut -d "." -f1-1 )"
    R2_name="$( echo $R2 | cut -d "." -f1-1 )"

    cutadapt -q 30,30 -u 20 -a $adapt1 -A $adapt2 -g file:"${folder}barcodes_unique.fa" \
              --minimum-length 20 -o "${R1_name}_clean.fastq.gz" \
              -p "${R2_name}_clean.fastq.gz" $R1 $R2

  elif [[ $folder == *"_2015_SoilGenomes"* ]]; then
    continue
    adapt1="AGATCGGAAGAGC"
    adapt2="AGATCGGAAGAGC"
    R1="${folder}"*"_R1.fastq"
    R2="${folder}"*"_R2.fastq"

    grep ^@HWI $R1 | cut -d":" -f10 | sort -r | uniq > "${folder}"barcodes_R1.txt
    grep ^@HWI $R2 | cut -d":" -f10 | sort -r | uniq > "${folder}"barcodes_R2.txt
    cat "${folder}barcodes_R"?".txt" >> "${folder}barcodes.txt"
    sort -u "${folder}barcodes.txt" >> "${folder}barcodes_unique.txt"
    python /N/dc2/projects/muri2/Task2/LTDE/bash/add_fasta_headers.py "${folder}barcodes_unique.txt" "${folder}barcodes_unique.fa"

    R1_name="$( echo $R1 | cut -d "." -f1-1 )"
    R2_name="$( echo $R2 | cut -d "." -f1-1 )"

    cutadapt -q 30,30 -u 20 -a $adapt1 -A $adapt2 -g file:"${folder}barcodes_unique.fa" \
              --minimum-length 20 -o "${R1_name}_clean.fastq" \
              -p "${R2_name}_clean.fastq" $R1 $R2

  elif [[ $folder == *"_GSF911"* ]]; then
    #continue
    R1="${folder}"*"_R1_001.fastq.gz"
    R2="${folder}"*"_R2_001.fastq.gz"
    # no barcode listed, remove overrepresented sequences after
    #zgrep ^@M01529 $R1 | cut -d":" -f10 | sort -r | uniq > "${folder}"barcodes_R1.txt
    #zgrep ^@M01529 $R2 | cut -d":" -f10 | sort -r | uniq > "${folder}"barcodes_R2.txt

    #M01529
    R1_name="$( echo $R1 | cut -d "." -f1-1 )"
    R2_name="$( echo $R2 | cut -d "." -f1-1 )"
    adapt=/N/dc2/projects/muri2/Task2/LTDE/bash/GSF911_adapters.fasta
    cutadapt -q 30,30 -u 20 -a file:"${adapt}" \
              --minimum-length 20 -o "${R1_name}_clean.fastq" \
              -p "${R2_name}_clean.fastq" $R1 $R2
  else
    continue

  #fastqc $R1 > "${folder}"R1_results.out
  #fastqc $R2 > "${folder}"R2_results.out

  #fastqc "${R1_name}_clean.fastq.gz" > "${folder}"R1_clean_results.out
  #fastqc "${R2_name}_clean.fastq.gz" > "${folder}"R2_clean_results.out
  fi
done


#fastqc Ev712_TAAGGCGAATCT-TATCCTCT_L002_R1_001_clean.fastq.gz > KBS0712_R1_clean_results.out


#zcat Ev13483_TAAGGCGAATCT-AAGGAGTA_L002_R1_001.fastq.gz | head -n 1
