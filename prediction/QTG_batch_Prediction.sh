#!/bin/env bash

#SBATCH -J causalQTL
#SBATCH -p DPB
#SBATCH -c 24
#SBATCH --mem=0

module purge;
module load Python/3.6.0 

srun python QTG_Finder.py Arabidopsis_features-v3.05.txt SSQ_batch_QTL_genes.csv

