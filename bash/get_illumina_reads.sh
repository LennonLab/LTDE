#!/bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=8,vmem=100gb,walltime=12:00:00
#PBS -M wrshoema@iu.edu
#PBS -m abe
#PBS -j oe

module unload python
module load python/3.6.1
module load fastqc
module load cutadapt

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

for folder in /N/dc2/projects/muri2/Task2/LTDE/illumina_data/*/
do
  echo $folder
done

for folder in /N/dc2/projects/muri2/Task2/LTDE/illumina_data/*/
do
  if [[ $folder == *"D400_100"* ]]; then
    sample_sheet="${folder}SampleSheet.csv"
    adapt1="$(sed '2q;d' ${sample_sheet} | cut -d"," -f5 | cut -d"-" -f1)"
    adapt2="$(sed '2q;d' ${sample_sheet} | cut -d"," -f5 | cut -d"-" -f2)"
    R1="${folder}"*"_R1_001.fastq.gz"
    R2="${folder}"*"_R2_001.fastq.gz"
    zgrep ^@HSQ $R1 | cut -d":" -f10 | sort -r | uniq > "${folder}"barcodes_R1.txt
    zgrep ^@HSQ $R2 | cut -d":" -f10 | sort -r | uniq > "${folder}"barcodes_R2.txt
    cat "${folder}barcodes_R"?".txt" >> "${folder}barcodes.txt"
    sort -u "${folder}barcodes.txt" >> "${folder}barcodes_unique.txt"
    python /N/dc2/projects/muri2/Task2/LTDE/bash/add_fasta_headers.py "${folder}barcodes_unique.txt" "${folder}barcodes_unique.fa"

    R1_name="$( echo "${R1}" | cut -d "." -f1-1 )"
    R2_name="$( echo "${R2}" | cut -d "." -f1-1 )"

    cutadapt -q 30,30 -a $adapt1 -A $adapt2 -g file:"${folder}barcodes_unique.fa" \
              --minimum-length 20 -o "${R1_name}_clean.fastq.gz" \
              -p "${R2_name}_clean.fastq.gz" $R1 $R2



  elif [[ $folder == *"_2015_SoilGenomes"* ]]; then
    adapt1="AGATCGGAAGAGC"
    adapt2="AGATCGGAAGAGC"
    R1="${folder}"*"_R1.fastq"
    R2="${folder}"*"_R2.fastq"
    # make sure this works
    grep ^@HSQ $R1 | cut -d":" -f10 | sort -r | uniq > "${folder}"barcodes_R1.txt
    grep ^@HSQ $R2 | cut -d":" -f10 | sort -r | uniq > "${folder}"barcodes_R2.txt
    cat "${folder}barcodes_R"?".txt" >> "${folder}barcodes.txt"
    sort -u "${folder}barcodes.txt" >> "${folder}barcodes_unique.txt"
    python /N/dc2/projects/muri2/Task2/LTDE/bash/add_fasta_headers.py "${folder}barcodes_unique.txt" "${folder}barcodes_unique.fa"

    R1_name="$( echo "${R1}" | cut -d "." -f1-1 )"
    R2_name="$( echo "${R2}" | cut -d "." -f1-1 )"

    cutadapt -q 30,30 -a $adapt1 -A $adapt2 -g file:"${folder}barcodes_unique.fa" \
              --minimum-length 20 -o "${R1_name}_clean.fastq" \
              -p "${R2_name}_clean.fastq" $R1 $R2

  elif [[ $folder == *"_GSF911"* ]]; then
    echo "What"
    # use GSF911_adapters.txt
    # cutadapt -a file:/N/dc2/projects/muri2/Task2/LTDE/bash/GSF911_adapters.fasta


  else
    echo "Unknown file"


  fastqc "$R1" --outdir="${folder}R1_fastqc"
  fastqc "$R2" --outdir="${folder}R2_fastqc"

  fastqc "${R1_name}_clean.fastq.gz" --outdir="${R1_name}_clean"
  fastqc "${R2_name}_clean.fastq.gz" --outdir="${R2_name}_clean"

  fi
done





#zcat Ev13483_TAAGGCGAATCT-AAGGAGTA_L002_R1_001.fastq.gz | head -n 1




# run quast after assembly





#fastqc $R1 >> results.out
#fastqc $R2 >> results.out

#cutadapt -q 10 -a AGATCGGAAGAGC --minimum-length 20 -o ./tmp.1.fastq -p ./tmp.2.fastq ./$R1 ./$R2
#cutadapt -q 10 -a AGATCGGAAGAGC --minimum-length 20 -o ./$GENOME.R2.trim.fastq -p ./$GENOME.R1.trim.fastq ./tmp.2.fastq ./tmp.1.fastq
#rm ./tmp.1.fastq ./tmp.2.fastq
#interleave-reads.py ./$GENOME.R1.trim.fastq ./$GENOME.R2.trim.fastq -o ./$GENOME.trim.interleaved.fastq
#fastq_quality_filter -Q33 -q 30 -p 50 -i ./$GENOME.trim.interleaved.fastq > ./$GENOME.trim.interleaved.trim.fastq
#extract-paired-reads.py ./$GENOME.trim.interleaved.trim.fastq
#normalize-by-median.py -k 25 -C 25 -N 4 -x 2e9 -p ./$GENOME.trim.interleaved.trim.fastq.pe
#extract-paired-reads.py ./$GENOME.trim.interleaved.trim.fastq.pe.keep
