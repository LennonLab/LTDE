#!/bin/bash

breseq_essentials="/N/dc2/projects/muri2/Task2/LTDE/data/breseq_essentials"
mkdir -p $breseq_essentials

for folder in "/N/dc2/projects/muri2/Task2/LTDE/data/breseq_output/"*
do
  strain="$(  echo "$folder" | cut -d"/" -f10-10)"
  mkdir -p "${breseq_essentials}/$strain"
  annotated="${folder}/output/evidence/annotated.gd"
  evidence="${folder}/output/evidence/evidence.gd"
  output="${folder}/output/output.gd"
  cp $annotated "${breseq_essentials}/$strain/annotated.gd"
  cp $evidence "${breseq_essentials}/$strain/evidence.gd"
  cp $output "${breseq_essentials}/$strain/output.gd"
done
