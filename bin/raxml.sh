#!/bin/bash
#PBS -k o
#PBS -l nodes=2:ppn=8,vmem=100gb,walltime=240:00:00
#PBS -M wrshoema@indiana.edu
#PBS -m abe
#PBS -j oe

# version 8.2.11
#module load raxml

raxmlHPC-PTHREADS -T 4 -f a -m GTRGAMMA -p 12345 -x 12345 -o NC_005042.1_353331-354795 -# autoMRE \
    -s /Users/WRShoemaker/GitHub/LTDE/data/align/ltde_seqs.good.filter.fasta \
    -n ltde_seqs -w /Users/WRShoemaker/GitHub/LTDE/data/tree

# for the protein tree..........
raxmlHPC-PTHREADS -T 4 -f a -m PROTGAMMALG -p 12345 -x 12345 -# autoMRE \
    -s \
    -n \
    -w

# -T = number of threads
# -f = specifies bootstrapping algorithm with ML generating tree at same time
# -m = substitution model, generalized time reversible gamma
# -p = starts tree randomly
# -x = starts tree randomly
# -o = outgroup (name after fasta entry)
# -#  = determines number of bootstrap replicates
# -s = aligned fasta file to be analyzed
# -n = name of output file
