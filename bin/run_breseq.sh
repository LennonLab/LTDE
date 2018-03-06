#!/bin/bash

OUT_breseq="/N/dc2/projects/muri2/Task2/LTDE/data/breseq_output/"
mkdir -p $OUT_breseq
breseq_bash="/N/dc2/projects/muri2/Task2/LTDE/bin/breseq/"
mkdir -p $breseq_bash

for folder in "/N/dc2/projects/muri2/Task2/LTDE/data/reads_clean/"*
do
  R1="${folder}/"*_R1_paired.fastq.gz
  R2="${folder}/"*_R2_paired.fastq.gz
  strain="$(  echo "$folder" | cut -d"/" -f10-10 | cut -d"-" -f2-2 )"
  strain_rep="$(  echo "$folder" | cut -d"/" -f10-10 | cut -d"-" -f2-3 )"
  if [ "${strain_rep}" = "ATCC43928-C1" ]; then
    C1="/N/dc2/projects/muri2/Task2/LTDE/data/reads_clean/GSF1046-ATCC43928-C2/GSF1046-ATCC43928-C2_clean_R1_paired.fastq.gz"
    C2="/N/dc2/projects/muri2/Task2/LTDE/data/reads_clean/GSF1046-ATCC43928-C2/GSF1046-ATCC43928-C2_clean_R2_paired.fastq.gz"
    reads="${R1} ${R2} ${C1} ${C2}"
    bash_out="${breseq_bash}/ATCC43928-C_breseq.sh"
    OUT_breseq_strain="${OUT_breseq}ATCC43928-C"
  elif [ "${strain_rep}" = "ATCC43928-C2" ]; then
    continue
  else
    reads="${R1} ${R2}"
    bash_out="${breseq_bash}${strain_rep}_breseq.sh"
    OUT_breseq_strain="${OUT_breseq}${strain_rep}"
  fi
  mkdir -p $OUT_breseq_strain
  ref="/N/dc2/projects/muri2/Task2/LTDE/data/reference_genomes/"*"/${strain}/G-Chr1.gbk"

  if [ -f $bash_out ]; then
    rm $bash_out
  fi

  echo '#!/bin/bash' >> $bash_out
  echo '#PBS -k o' >> $bash_out
  echo '#PBS -l nodes=1:ppn=8,vmem=100gb,walltime=30:00:00' >> $bash_out
  echo '#PBS -M wrshoema@umail.iu.edu' >> $bash_out
  echo '#PBS -m abe' >> $bash_out
  echo '#PBS -j oe' >> $bash_out
  echo '' >> $bash_out
  echo 'module load breseq' >> $bash_out
  echo "breseq -j 8 -p -o ${OUT_breseq_strain} -r ${ref} ${reads}" >> $bash_out
  qsub $bash_out
done

# code to get reference genomes...
