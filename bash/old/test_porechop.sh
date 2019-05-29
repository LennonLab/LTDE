#!/bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=8,vmem=100gb,walltime=16:00:00
#PBS -M wrshoema@iu.edu
#PBS -m abe
#PBS -j oe

module unload python
module load python/3.6.1

PC=/N/u/wrshoema/Carbonate/Porechop/porechop-runner.py
IN2=/N/dc2/projects/muri2/Task2/LTDE/data/nanopore_basecalled/run2
IN3=/N/dc2/projects/muri2/Task2/LTDE/data/nanopore_basecalled/run3
IN4=/N/dc2/projects/muri2/Task2/LTDE/data/nanopore_basecalled/run4

BARCODE_DIR2=/N/dc2/projects/muri2/Task2/LTDE/data/nanopore_basecalled_pcbc/run2
BARCODE_DIR3=/N/dc2/projects/muri2/Task2/LTDE/data/nanopore_basecalled_pcbc/run3
BARCODE_DIR4=/N/dc2/projects/muri2/Task2/LTDE/data/nanopore_basecalled_pcbc/run4

#$PC -i $IN2 --format fastq \
#      -b $BARCODE_DIR2 \
#      --barcode_threshold 75 \
#      --barcode_diff 5 \
#      --require_two_barcodes \
#      --discard_middle

#$PC -i $IN3 --format fastq \
#      -b $BARCODE_DIR3 \
#      --barcode_threshold 75 \
#      --barcode_diff 5 \
#      --require_two_barcodes \
#      --discard_middle

$PC -i $IN4 --format fastq \
      -b $BARCODE_DIR4 \
      --barcode_threshold 75 \
      --barcode_diff 5 \
      --require_two_barcodes \
      --discard_middle
