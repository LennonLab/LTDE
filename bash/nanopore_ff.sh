#!/bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=8,vmem=100gb,walltime=320:00:00
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

cfg=/N/dc2/projects/muri2/Task2/LTDE/ont-guppy-cpu/data/dna_r9.4.1_450bps_flipflop.cfg
fast5_dir_run2=/N/dc2/projects/muri2/Task2/LTDE/GSF2099/GSF2099-Lennon-run2-20190326/DNA-NativeBarcode/GSF2099-set2/20190326_2329_MN22062_FAK62803_d8b2ba25/fast5_pass
out_dir_run2=/N/dc2/projects/muri2/Task2/LTDE/data/nanopore_basecalled_ff/run2

/N/dc2/projects/muri2/Task2/LTDE/ont-guppy-cpu/bin/guppy_basecaller \
    -i $fast5_dir_run2 \
    -t 12 \
    -s $out_dir_run2 \
    -k SQK-LSK109 \
    -f FLO-MIN106 \
    --config $cfg


fast5_dir_run3=/N/dc2/projects/muri2/Task2/LTDE/GSF2099/GSF2099-Lennon-run3-20190408/GSF2099-set3/20190408_2058_MN22062_FAK70247_d34ce651/fast5_pass
out_dir_run3=/N/dc2/projects/muri2/Task2/LTDE/data/nanopore_basecalled_ff/run3
/N/dc2/projects/muri2/Task2/LTDE/ont-guppy-cpu/bin/guppy_basecaller \
    -i $fast5_dir_run3 \
    -t 12 \
    -s $out_dir_run3 \
    -k SQK-LSK109 \
    -f FLO-MIN106 \
    --config $cfg

fast5_dir_run4=/N/dc2/projects/muri2/Task2/LTDE/GSF2099-set4/1/20190422_1815_MN22062_FAK67864_08d7efdf/fast5_pass
out_dir_run4=/N/dc2/projects/muri2/Task2/LTDE/data/nanopore_basecalled_ff/run4
/N/dc2/projects/muri2/Task2/LTDE/ont-guppy-cpu/bin/guppy_basecaller \
    -i $fast5_dir_run4 \
    -t 12 \
    -s $out_dir_run4 \
    -k SQK-LSK109 \
    -f FLO-MIN106 \
    --config $cfg
