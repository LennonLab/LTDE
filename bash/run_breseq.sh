#!/bin/bash

# no 703 for now
#"ATCC13985-A"
declare -a samples=("ATCC13985-B" "ATCC13985-C" "ATCC13985-D" \
                    "ATCC43928-A" "ATCC43928-B" "ATCC43928-C" "ATCC43928-D" \
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


OUT_breseq="/N/dc2/projects/muri2/Task2/LTDE/data/breseq/"
mkdir -p $OUT_breseq
breseq_bash="/N/dc2/projects/muri2/Task2/LTDE/bin/breseq/"
mkdir -p $breseq_bash
mkdir -p /N/dc2/projects/muri2/Task2/LTDE/data/breseq_out/
mkdir -p /N/dc2/projects/muri2/Task2/LTDE/data/breseq_err/


for sample in "${samples[@]}"
do
  taxon="$(  echo "$sample" | cut -d"-" -f1-1 )"
  reads="/N/dc2/projects/muri2/Task2/LTDE/data/reads_clean/"*"${sample}"*
  bash_script="/N/dc2/projects/muri2/Task2/LTDE/bash/breseq_scripts/${sample}.sh"
  if [ -f $bash_script ]; then
    rm $bash_script
  fi

  OUT_breseq_out="/N/dc2/projects/muri2/Task2/LTDE/data/breseq_out/${sample}.out"
  OUT_breseq_err="/N/dc2/projects/muri2/Task2/LTDE/data/breseq_err/${sample}.err"
  OUT_breseq="/N/dc2/projects/muri2/Task2/LTDE/data/breseq/${sample}"

  ref="/N/dc2/projects/muri2/Task2/LTDE/data/genomes_ncbi/${taxon}/"*".gbff"

  echo '#!/bin/bash' >> $bash_script
  echo '#PBS -k o' >> $bash_script
  echo '#PBS -l nodes=1:ppn=8,vmem=100gb,walltime=24:00:00' >> $bash_script
  echo '#PBS -M wrshoema@iu.edu' >> $bash_script
  echo '#PBS -m abe' >> $bash_script
  echo '#PBS -j oe' >> $bash_script
  echo '' >> $bash_script
  echo 'module load breseq' >> $bash_script
  echo "breseq -j 8 -p -o ${OUT_breseq} -r ${ref} ${reads} > ${OUT_breseq_out} 2> ${OUT_breseq_err}" >> $bash_script
  qsub $bash_script
done
