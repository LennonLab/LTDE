#!/bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=8,vmem=100gb,walltime=1:00:00
#PBS -M wrshoema@indiana.edu
#PBS -m abe
#PBS -j oe

#module load gcc/4.9.2
#module load boost/1.52.0
#module load openmpi
#module load mothur/1.38.1
module load mothur


cd /N/dc2/projects/Lennon_Sequences/LTDE_Tree

#wget https://www.mothur.org/w/images/d/dc/Trainset16_022016.rdp.tgz
#wget https://www.mothur.org/w/images/7/71/Silva.seed_v132.tgz
#tar -xzvf Silva.seed_v132.tgz
#mv Trainset16_022016.rdp.tgz silva.v4.fasta



mothur persistence_707_712_721.batch
