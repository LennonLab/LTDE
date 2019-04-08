#!/bin/bash
#PBS -k o
#PBS -l nodes=2:ppn=8,vmem=100gb,walltime=240:00:00
#PBS -M wrshoema@indiana.edu
#PBS -m abe
#PBS -j oe

module load epa

cd /N/dc2/projects/Lennon_Sequences/LTDE_Tree

epa-ng -T 4 -s ./nmicrobiol201648-s7_clean.txt -t ./nmicrobiol201648-s8_clean.txt -q ./persistence_707_712_721_test.fasta -w ./epa_ng
