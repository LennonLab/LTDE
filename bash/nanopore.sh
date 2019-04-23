#!/bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=8,vmem=100gb,walltime=72:00:00
#PBS -M wrshoema@iu.edu
#PBS -m abe
#PBS -j oe

# base-calling with Guppy

# FLO-MIN106
# compatible flowcells
#Pore version: R9.4.1 FLO-MIN106
#Barcoding kit ID: 1D Native barcoding genomic DNA, with EXP-NBD104 barcode
#Library ligation kit ID: SQK-LSK109
# dna_r9.4.1_450bps
#<strand type>_<pore version>_<speed>[custom tags].cfg

# --calib_detect
# -k name of kit
# -f name of flowcell

cfg=/N/dc2/projects/muri2/Task2/LTDE/ont-guppy-cpu/data/dna_r9.4.1_450bps.cfg
fast5_dir_run2=/N/dc2/projects/muri2/Task2/LTDE/GSF2099/GSF2099-Lennon-run2-20190326/DNA-NativeBarcode/GSF2099-set2/20190326_2329_MN22062_FAK62803_d8b2ba25/fast5_pass
out_dir_run2=/N/dc2/projects/muri2/Task2/LTDE/data/nanopore_basecalled/run2

fast5_dir_run3=/N/dc2/projects/muri2/Task2/LTDE/GSF2099/GSF2099-Lennon-run3-20190408/GSF2099-set3/20190408_2058_MN22062_FAK70247_d34ce651/fast5_pass
out_dir_run3=/N/dc2/projects/muri2/Task2/LTDE/data/nanopore_basecalled/run3

#fast5_dir_run2=/N/dc2/projects/muri2/Task2/LTDE/GSF2099/test

#/N/dc2/projects/muri2/Task2/LTDE/ont-guppy-cpu/bin/guppy_basecaller \
#    -i $fast5_dir_run2 \
#    -t 12 \
#    -s $out_dir_run2 \
#    -k SQK-LSK109 \
#    -f FLO-MIN106 \
#    --config $cfg

#/N/dc2/projects/muri2/Task2/LTDE/ont-guppy-cpu/bin/guppy_basecaller \
#    -i $fast5_dir_run3 \
#    -t 12 \
#    -s $out_dir_run3 \
#    -k SQK-LSK109 \
#    -f FLO-MIN106 \
#    --config $cfg



#Barcoding/demultiplexing within specific kits
out_dir_bc_run2=/N/dc2/projects/muri2/Task2/LTDE/data/nanopore_basecalled_bc/run2
out_dir_bc_run3=/N/dc2/projects/muri2/Task2/LTDE/data/nanopore_basecalled_bc/run3

bc_cfg=/N/dc2/projects/muri2/Task2/LTDE/ont-guppy-cpu/data/barcoding/configuration.cfg
/N/dc2/projects/muri2/Task2/LTDE/ont-guppy-cpu/bin/guppy_barcoder \
    --input_path $out_dir_run2
    --save_path $out_dir_bc_run2
    --config $bc_cfg


#/N/dc2/projects/muri2/Task2/LTDE/ont-guppy-cpu/bin/guppy_barcoder \
#    --input_path $out_dir_run3
#    --save_path $out_dir_bc_run3
#    --config $bc_cfg

# clean fasta reads
# then assemble with spades
