#!/bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=8,vmem=100gb,walltime=1000:00:00
#PBS -M wrshoema@indiana.edu
#PBS -m abe
#PBS -j oe

# version 8.2.11
module load raxml

#raxmlHPC-PTHREADS -T 4 -f a -m GTRGAMMA -p 12345 -x 12345 -o NC_005042.1_353331-354795 -# autoMRE \
#    -s /Users/WRShoemaker/GitHub/LTDE/data/align/16S_rRNA_seqs/ltde_seqs.good.filter.fasta \
#    -n ltde_seqs_test -w /Users/WRShoemaker/GitHub/LTDE/data/tree

# test for the protein tree..........
raxmlHPC-PTHREADS -T 4 -f a -m PROTGAMMALG -p 12345 -x 12345 -# autoMRE \
    -s /Users/WRShoemaker/GitHub/LTDE/data/align/ribosomal_protein_seqs_align_concat.fa \
    -n test_protein \
    -w /Users/WRShoemaker/GitHub/LTDE/data/align

raxmlHPC-PTHREADS -T 4 -f a -m PROTGAMMALG -p 12345 -x 12345 -# autoMRE \
    -s /N/dc2/projects/muri2/Task2/LTDE/data/align/ribosomal_protein_seqs_align_concat.fa \
    -n protein_tree \
    -w /N/dc2/projects/muri2/Task2/LTDE/data/tree


# -T = number of threads
# -f = specifies bootstrapping algorithm with ML generating tree at same time
# -m = substitution model, generalized time reversible gamma
# -p = starts tree randomly
# -x = starts tree randomly
# -o = outgroup (name after fasta entry)
# -#  = determines number of bootstrap replicates
# -s = aligned fasta file to be analyzed
# -n = name of output file
