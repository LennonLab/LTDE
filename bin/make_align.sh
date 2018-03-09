#!/bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=8,vmem=100gb,walltime=48:00:00
#PBS -M wrshoema@indiana.edu
#PBS -m abe
#PBS -j oe


module load mothur


cd /N/dc2/projects/Lennon_Sequences/LTDE_Tree

#wget https://www.mothur.org/w/images/d/dc/Trainset16_022016.rdp.tgz
#wget https://www.mothur.org/w/images/7/71/Silva.seed_v132.tgz
#tar -xzvf Silva.seed_v132.tgz
#mv Trainset16_022016.rdp.tgz silva.v4.fasta



mothur make_align.batch
