#!/usr/bin/python

# Version 2.0
# last update: 20200131
# Purpose: This code is used for re-train the QTG models if users want to add features or use different parameters for modeling. 
# Usage = "QTG_Finder_train.py -fl feature list -sp species_abbreviation"
# feature list: Arabidopsis_features_v5_ortholog.csv for Arabidopsis; rice_features_v5_ortholog.csv for rice ; Setaria_features.csv for Setaria viridis; Sorghum_features.csv for Sorghum 
# species_abbreviation: "AT" for Arabidopsis; "OS" for rice; "SB" for sorghum; "SV" for Setaria
# Usage example: python3 QTG_Finder_train.py -fl Arabidopsis_features_v5_ortholog.csv -sp 'AT'
# Usage example: python3 QTG_Finder_train.py -fl Setaria_features.csv -sp 'SV' --rmgene Sevir.9G121800,Sevir.5G410400

import numpy as np
import random
import pandas as pd
from sklearn import preprocessing
from sklearn import ensemble
import time
import argparse
import pickle

parser = argparse.ArgumentParser(description='QTL_Finder_v2.0')
requiredArg = parser.add_argument_group('required arguments')
requiredArg.add_argument('-fl',dest='fl',required=True,help='Input feature list: use Arabidopsis_features_v4.csv for Arabidopsis; use rice_features_v2.csv for rice')
requiredArg.add_argument('-sp',dest='sp',required=True,choices=['AT','OS','SV','SB'], help='Species: use "AT" for Arabidopsis; "OS" for rice; "SV" for Setaria viridis; "SB" for Sorghum bicolor')
parser.add_argument("--rmgene",dest='gene_removal',help="remove specific genes from training set,seperate genes by comma, only used for test")
args = parser.parse_args()

def train_qtg(df, train_set):
    random.seed(11)
    if args.sp=='AT': 
        mol_para=['auto',20,2]# parameters for Arabidopsis
    if args.sp=='OS':
        mol_para=['auto',5,2]#  parameters for rice
    if args.sp=='SV':
        mol_para=['auto',5,6] # parameter for Setaria
    if args.sp=='SB':
        mol_para=['auto',5,6] # parameter for Sorghum
    clf = ensemble.RandomForestClassifier(n_estimators=200, min_samples_split=mol_para[2],max_features=mol_para[0],n_jobs=-1) # Random forest parameters
    neg_inter=5000  # interrations for randomly selecting negatives from genome and re-training models 
    pickle.dump(neg_inter, pik_f)
    for i in range(0, neg_inter):
                train_data = train_set # positives used for training 
                training_negative = random.sample(list(df[df['class']==0].index), int(len(train_data)*mol_para[1])) # randomly select negatives from genome genes 
                train_data=train_data.append(df.iloc[training_negative]) 
                train_feature=train_data.drop(['class'], axis=1) 
                clf.fit(train_feature, train_data['class'])  # model fitting
                pickle.dump(clf, pik_f)
    print('training complete')
    
#### check known QTL gene on QTLs, remove, reducant ID, run
if __name__ == '__main__':
    with open("../"+args.sp+"_model.dat", "wb") as pik_f:
        dt = args.fl # input feature list
        start_time = time.time()
        df = pd.read_csv(dt)
        if args.gene_removal!=None:
            rm_gene_entry=args.gene_removal
            rm_gene_list=rm_gene_entry.split(",")
            if len(rm_gene_list)>0: 
                for i in rm_gene_list: 
                    index=df[df['ID']==i].index[0]
                    if df.iloc[index,-1]==1: 
                        df.iloc[index,-1]=0
        pickle.dump(df, pik_f)
        df=df.dropna(axis=1,how='all') 
        df = df.drop(['ID'], axis=1) 
        train_set = df[df['class']==1]
        train_qtg(df, train_set) # execute 

    print("--- %s seconds ---" % round((time.time() - start_time),2) ) 

