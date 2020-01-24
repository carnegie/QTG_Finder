#!/bin/env bash

#SBATCH -J sorgh_causalQTL
#SBATCH -p DPB
#SBATCH -c 24
#SBATCH --mem=0

module purge;
module load Python/3.6.0 

srun python QTG_full_model_testing_FAN_1_1.py sobic_fea5.csv ./input/SB/PRR37.csv 'SB'
srun python QTG_full_model_testing_FAN_1_1.py sobic_fea5.csv ./input/SB/PhyB_ma3_adj_cord_causal_cent_substraction.csv 'SB'
srun python QTG_full_model_testing_FAN_1_1.py sobic_fea5.csv ./input/SB/Ghd7_Ma6.csv 'SB'
srun python QTG_full_model_testing_FAN_1_1.py sobic_fea5.csv ./input/SB/sbMATE1_.csv 'SB'
srun python QTG_full_model_testing_FAN_1_1.py sobic_fea5.csv ./input/SB/ds1.csv 'SB'
srun python QTG_full_model_testing_FAN_1_1.py sobic_fea5.csv ./input/SB/Sh1.csv 'SB'
srun python QTG_full_model_testing_FAN_1_1.py sobic_fea5.csv ./input/SB/wx.csv 'SB'
srun python QTG_full_model_testing_FAN_1_1.py sobic_fea5.csv ./input/SB/dw2.csv 'SB'
srun python QTG_full_model_testing_FAN_1_1.py sobic_fea5.csv ./input/SB/rf2_PPR.csv 'SB'
srun python QTG_full_model_testing_FAN_1_1.py sobic_fea5.csv ./input/SB/bmr2.csv 'SB'

