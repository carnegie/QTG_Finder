#!/bin/env bash

#SBATCH -J causalQTL
#SBATCH -p DPB
#SBATCH -c 24
#SBATCH --mem=0

module purge;
module load Python/3.6.0 

srun python QTG_full_model_testing_FAN_1_1.py Arabidopsis_features_v5_only_ortho.csv ./input/AT/Huang_2012.csv 'AT'
srun python QTG_full_model_testing_FAN_1_1.py Arabidopsis_features_v5_only_ortho.csv ./input/AT/Fbien_QTL_genes.csv 'AT'
srun python QTG_full_model_testing_FAN_1_1.py Arabidopsis_features_v5_only_ortho.csv ./input/AT/Motte_2014.csv 'AT'
srun python QTG_full_model_testing_FAN_1_1.py Arabidopsis_features_v5_only_ortho.csv ./input/AT/yuan_2016_XQTL_list.csv 'AT'
srun python QTG_full_model_testing_FAN_1_1.py Arabidopsis_features_v5_only_ortho.csv ./input/AT/Guo_2016_SSQ3_AHK3.csv 'AT'
srun python QTG_full_model_testing_FAN_1_1.py Arabidopsis_features_v5_only_ortho.csv ./input/AT/Guo_2016_SSQ10_AHK2.csv 'AT'
srun python QTG_full_model_testing_FAN_1_1.py Arabidopsis_features_v5_only_ortho.csv ./input/AT/TGG.csv 'AT'
srun python QTG_full_model_testing_FAN_1_1.py Arabidopsis_features_v5_only_ortho.csv ./input/AT/GSOH1.csv 'AT'


