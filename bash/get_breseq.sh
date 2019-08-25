#!/bin/bash

breseq_essentials="/N/dc2/projects/muri2/Task2/LTDE/data/breseq_essentials"
mkdir -p /N/dc2/projects/muri2/Task2/LTDE/data/breseq/evidence
mkdir -p /N/dc2/projects/muri2/Task2/LTDE/data/breseq/annotated
mkdir -p /N/dc2/projects/muri2/Task2/LTDE/data/breseq/output
mkdir -p /N/dc2/projects/muri2/Task2/LTDE/data/breseq/summary

for folder in "/N/dc2/projects/muri2/Task2/LTDE/data/breseq/"*
do
  strain="$(  echo "$folder" | cut -d"/" -f10-10)"
  #cp "${folder}/output/evidence/annotated.gd" "/N/dc2/projects/muri2/Task2/LTDE/data/breseq/annotated/${strain}.gd"
  #cp "${folder}/output/evidence/evidence.gd" "/N/dc2/projects/muri2/Task2/LTDE/data/breseq/evidence/${strain}.gd"
  #cp "${folder}/output/output.gd" "/N/dc2/projects/muri2/Task2/LTDE/data/breseq/output/${strain}.gd"
  cp "${folder}/data/summary.json" "/N/dc2/projects/muri2/Task2/LTDE/data/breseq/summary/${strain}.json"
done
