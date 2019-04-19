#!/bin/bash

# base-calling with Guppy

# FLO-MIN106
# compatible flowcells
# dna_r9.4.1_450bps.cfg,
cfg=/N/dc2/projects/muri2/Task2/LTDE/ont-guppy-cpu/data/
#/N/dc2/projects/muri2/Task2/LTDE/ont-guppy-cpu/bin/guppy_basecaller --config dna_r9.4_450bps.cfg

--config $cfg -i $fast5_dir -t 12 -s $out_dir --enable_trimming
# --calib_detect
# -k name of kit
# -f name of flowcell

# dna_r9.4.1_450bps
#<strand type>_<pore version>_<speed>[custom tags].cfg

#Barcoding/demultiplexing within specific kits

# guppy_barcoder --input_path <folder containing fastq files>
# --save_path <output folder>
# --config /N/dc2/projects/muri2/Task2/LTDE/ont-guppy-cpu/data/barcoding/configuration.cfg



# clean fasta reads
# then assemble with spades
