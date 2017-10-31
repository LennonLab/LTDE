#!/bin/bash
#PBS -k o
#PBS -l nodes=2:ppn=8,vmem=100gb,walltime=5:00:00
#PBS -M wrshoema@indiana.edu
#PBS -m abe
#PBS -j oe

module load raxml/8.0.26

# cd into the directory with your alignment

cd /N/dc2/projects/Lennon_Sequences/LTDE_Tree/tree

raxmlHPC-PTHREADS -T 4 -f a -m GTRGAMMA -p 12345 -x 12345 -o Methanosarcina -# autoMRE -s ./persistence_707_712_721.good.filter.fasta -n LTDE.ML

# -T = number of threads
# -f = specifies bootstrapping algorithm with ML generating tree at same time
# -m = substitution model, generalized time reversible gamma
# -p = starts tree randomly
# -x = starts tree randomly
# -o = outgroup (name after fasta entry)
# -#  = determines number of bootstrap replicates
# -s = aligned fasta file to be analyzed
# -n = name of output file 