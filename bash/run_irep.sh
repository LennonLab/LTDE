#!/bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=8,vmem=100gb,walltime=72:00:00
#PBS -M wrshoema@indiana.edu
#PBS -m abe
#PBS -j oe

# iRep and bPTR uses paired-end information, which breseq doesn't use
# so we'll need to map the reads using bwa and then clean them a little
# with samtools

module load bwa
module load samtools
module unload python
module load python/3.6.1

# no 703 for now
#"ATCC13985-A"
#"ATCC43928-C"
declare -a samples=("ATCC13985-B" "ATCC13985-C" "ATCC13985-D" \
                    "ATCC43928-A" "ATCC43928-B" "ATCC43928-D" \
                    "KBS0702-A" "KBS0702-B" "KBS0702-C" "KBS0702-D" "KBS0702-E" "KBS0702-F" \
                    "KBS0705-A" "KBS0705-B" "KBS0705-C" "KBS0705-D" \
                    "KBS0706-A" "KBS0706-B" "KBS0706-C" "KBS0706-D" \
                    "KBS0707-A" "KBS0707-B" "KBS0707-C" "KBS0707-D" \
                    "KBS0710-A" "KBS0710-B" "KBS0710-C" "KBS0710-D" \
                    "KBS0711-A" "KBS0711-C" "KBS0711-D" "KBS0711-K" \
                    "KBS0712-A" "KBS0712-B" "KBS0712-C" "KBS0712-D" \
                    "KBS0713-A" "KBS0713-B" "KBS0713-C" \
                    "KBS0715-A" "KBS0715-B" "KBS0715-C" "KBS0715-D" \
                    "KBS0721-A" "KBS0721-B" "KBS0721-C" "KBS0721-D" \
                    "KBS0722-A" "KBS0722-B" "KBS0722-C" "KBS0722-D" \
                    "KBS0724-A" "KBS0724-B" "KBS0724-C" "KBS0724-D" \
                    "KBS0727-A" "KBS0727-B" "KBS0727-C" "KBS0727-D" \
                    "KBS0801-A" "KBS0801-B" "KBS0801-C" "KBS0801-D" \
                    "KBS0802-A" "KBS0802-B" "KBS0802-C" "KBS0802-D" \
                    "KBS0812-A" "KBS0812-B" "KBS0812-C" "KBS0812-D")


mkdir -p /N/dc2/projects/muri2/Task2/LTDE/data/bwa_bam
mkdir -p /N/dc2/projects/muri2/Task2/LTDE/data/iRep


for sample in "${samples[@]}"
do
  taxon="$(  echo "$sample" | cut -d"-" -f1-1 )"
  ref="/N/dc2/projects/muri2/Task2/LTDE/data/genomes_ncbi/${taxon}/"*".fna"
  bwa index $ref
  samtools faidx $ref
  if [ $sample = "ATCC43928-C" ]; then
    R1_1="/N/dc2/projects/muri2/Task2/LTDE/data/reads_clean/GSF1046-ATCC43928-C1_S43_R1_001_clean.fastq.gz"
    R2_1="/N/dc2/projects/muri2/Task2/LTDE/data/reads_clean/GSF1046-ATCC43928-C1_S43_R2_001_clean.fastq.gz"
    R1_2="/N/dc2/projects/muri2/Task2/LTDE/data/reads_clean/GSF1046-ATCC43928-C2_S50_R1_001_clean.fastq.gz"
    R2_2="/N/dc2/projects/muri2/Task2/LTDE/data/reads_clean/GSF1046-ATCC43928-C2_S50_R2_001_clean.fastq.gz"
    bwa mem -t 4 $ref $R1_1 $R2_1 | samtools view -F 4 -bT $ref - \
        | samtools sort -o "/N/dc2/projects/muri2/Task2/LTDE/data/bwa_bam/ATCC43928-C1.bam"
    bwa mem -t 4 $ref $R1_2 $R2_2 | samtools view -F 4 -bT $ref - \
        | samtools sort -o "/N/dc2/projects/muri2/Task2/LTDE/data/bwa_bam/ATCC43928-C2.bam"

    #samtools merge "/N/dc2/projects/muri2/Task2/LTDE/data/bwa_bam/ATCC43928-C.bam" \
    #    "/N/dc2/projects/muri2/Task2/LTDE/data/bwa_bam/ATCC43928-C1.bam" \
    #    "/N/dc2/projects/muri2/Task2/LTDE/data/bwa_bam/ATCC43928-C2.bam"

    samtools merge - \
        "/N/dc2/projects/muri2/Task2/LTDE/data/bwa_bam/ATCC43928-C1.bam" \
        "/N/dc2/projects/muri2/Task2/LTDE/data/bwa_bam/ATCC43928-C2.bam" \
        | samtools sort - \
        | samtools view -h -o "/N/dc2/projects/muri2/Task2/LTDE/data/bwa_bam/ATCC43928-C.sam"

  else
    R1="/N/dc2/projects/muri2/Task2/LTDE/data/reads_clean/"*"${sample}"*"_R1_"*
    R2="/N/dc2/projects/muri2/Task2/LTDE/data/reads_clean/"*"${sample}"*"_R2_"*
    bwa mem -t 4 $ref $R1 $R2 | samtools view -F 4 -bT $ref - \
        | samtools sort - \
        | samtools view -h -o "/N/dc2/projects/muri2/Task2/LTDE/data/bwa_bam/${sample}.sam"

    #-o "/N/dc2/projects/muri2/Task2/LTDE/data/bwa_bam/${sample}.bam"

  fi

  #samtools sort - "/N/dc2/projects/muri2/Task2/LTDE/data/bwa_bam/${sample}.bam" \
  #    | samtools view -h -o "/N/dc2/projects/muri2/Task2/LTDE/data/bwa_bam/${sample}.sam"

  /N/u/wrshoema/Carbonate/.local/bin/iRep -f $ref -ff \
      -s "/N/dc2/projects/muri2/Task2/LTDE/data/bwa_bam/${sample}.sam" \
      -o "/N/dc2/projects/muri2/Task2/LTDE/data/iRep/${sample}"
done
