#!/bin/env bash

#SBATCH -J causalQTL
#SBATCH -p DPB
#SBATCH -c 24
#SBATCH --mem=0

module purge;
module load Python/3.6.0 


srun python QTG_literature_validation.py ./input/rice_features_v1.3.11.txt ./input/OS/Qi_2008.csv 'OS'
srun python QTG_literature_validation.py ./input/rice_features_v1.3.11.txt ./input/OS/Liu_2017.csv 'OS'
srun python QTG_literature_validation.py ./input/rice_features_v1.3.11.txt ./input/OS/Hu_2008.csv 'OS'
srun python QTG_literature_validation.py ./input/rice_features_v1.3.11.txt ./input/OS/Gao_2016.csv 'OS'
srun python QTG_literature_validation.py ./input/rice_features_v1.3.11.txt ./input/OS/Oikawa_2015.csv 'OS'
srun python QTG_literature_validation.py ./input/rice_features_v1.3.11.txt ./input/OS/Fan_2016.csv 'OS'
srun python QTG_literature_validation.py ./input/rice_features_v1.3.11.txt ./input/OS/Fukuoka_2014.csv 'OS'
srun python QTG_literature_validation.py ./input/rice_features_v1.3.11.txt ./input/OS/Fjellstrom_2004_pikh.csv 'OS'
srun python QTG_literature_validation.py ./input/rice_features_v1.3.11.txt ./input/OS/zeng_2013.csv 'OS'