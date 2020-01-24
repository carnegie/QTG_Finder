#!/bin/env bash

#SBATCH -J causalQTL
#SBATCH -p DPB
#SBATCH -c 8
#SBATCH --mem=8000

module purge;
module load Python/3.6.0 

srun python QTG_literature_validation.py ./input/Arabidopsis_features_v4.csv./input/AT/Huang_2012.csv 'AT'
srun python QTG_literature_validation.py ./input/Arabidopsis_features_v4.csv./input/AT/Fbien_QTL_genes.csv 'AT'
srun python QTG_literature_validation.py ./input/Arabidopsis_features_v4.csv./input/AT/Motte_2014.csv 'AT'
srun python QTG_literature_validation.py ./input/Arabidopsis_features_v4.csv./input/AT/yuan_2016_XQTL_list.csv 'AT'
srun python QTG_literature_validation.py ./input/Arabidopsis_features_v4.csv./input/AT/Guo_2016_SSQ3_AHK3.csv 'AT'
srun python QTG_literature_validation.py ./input/Arabidopsis_features_v4.csv./input/AT/Guo_2016_SSQ10_AHK2.csv 'AT'
srun python QTG_literature_validation.py ./input/Arabidopsis_features_v4.csv./input/AT/TGG.csv 'AT'
srun python QTG_literature_validation.py ./input/Arabidopsis_features_v4.csv./input/AT/GSOH1.csv 'AT'