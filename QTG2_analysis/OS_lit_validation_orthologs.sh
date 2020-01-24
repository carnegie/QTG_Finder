#!/bin/env bash

#SBATCH -J causalQTL
#SBATCH -p DPB
#SBATCH -c 24
#SBATCH --mem=0

module purge;
module load Python/3.6.0 


srun python QTG_full_model_testing_FAN_1_1.py rice_features_v5_only_ortholog.csv ./input/OS/Qi_2008.csv 'OS'
srun python QTG_full_model_testing_FAN_1_1.py rice_features_v5_only_ortholog.csv ./input/OS/Liu_2017.csv 'OS'
srun python QTG_full_model_testing_FAN_1_1.py rice_features_v5_only_ortholog.csv ./input/OS/Hu_2008.csv 'OS'
srun python QTG_full_model_testing_FAN_1_1.py rice_features_v5_only_ortholog.csv ./input/OS/Gao_2016.csv 'OS'
srun python QTG_full_model_testing_FAN_1_1.py rice_features_v5_only_ortholog.csv ./input/OS/Oikawa_2015.csv 'OS'
srun python QTG_full_model_testing_FAN_1_1.py rice_features_v5_only_ortholog.csv ./input/OS/Fan_2016.csv 'OS'
srun python QTG_full_model_testing_FAN_1_1.py rice_features_v5_only_ortholog.csv ./input/OS/Fukuoka_2014.csv 'OS'
srun python QTG_full_model_testing_FAN_1_1.py rice_features_v5_only_ortholog.csv ./input/OS/Fjellstrom_2004_pikh.csv 'OS'
srun python QTG_full_model_testing_FAN_1_1.py rice_features_v5_only_ortholog.csv ./input/OS/zeng_2013.csv 'OS'
srun python QTG_full_model_testing_FAN_1_1.py rice_features_v5_only_ortholog.csv ./input/OS/Dixit_2018.csv 'OS'