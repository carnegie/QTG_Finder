#!/usr/bin/python

# Version 1.1
# last update: 20190307
# Purpose: This code is used for re-train the QTG models if user want to add features or use different parameters for models. 
# Usage = "QTG_Finder_train.py -fl feature list -sp species_abbreviation"
# feature list: use Arabidopsis_features_v4.csv for Arabidopsis; use rice_features_v2.csv for rice 
# species_abbreviation: "AT" for Arabidopsis; "OS" for rice
# Usage example: python QTG_Finder_train.py -fl Arabidopsis_features_v4.csv -sp 'AT'

import numpy as np
import random
import pandas as pd
from sklearn import preprocessing
from sklearn import ensemble
import time
import argparse
import pickle

parser = argparse.ArgumentParser(description='QTL_Finder_v1.0')
requiredArg = parser.add_argument_group('required arguments')
requiredArg.add_argument('-fl',dest='fl',required=True,help='Input feature list: use Arabidopsis_features_v4.csv for Arabidopsis; use rice_features_v2.csv for rice')
requiredArg.add_argument('-sp',dest='sp',required=True,choices=['AT', 'OS'], help='Species: use "AT" for Arabidopsis; "OS" for rice')
args = parser.parse_args()

def train_qtg(df, train_set):
    random.seed(11)
    if args.sp=='AT': 
        mol_para=[9,20]# optimum parameters for Arabidopsis
    if args.sp=='OS':
        mol_para=[9,5]# optimum parameters for rice
    clf = ensemble.RandomForestClassifier(n_estimators=100, min_samples_split=2,max_features=mol_para[0],n_jobs=-1) # Random forest parameters
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
    with open(args.sp+"_model.dat", "wb") as pik_f:
        dt = args.fl # input feature list
        start_time = time.time()
        df = pd.read_csv(dt)
        pickle.dump(df, pik_f)
        df=df.dropna(axis=1,how='all') 
        df['network_weight'] = preprocessing.scale(df['network_weight']) # standardization
        df = df.drop(['ID'], axis=1) 
        train_set = df[df['class']==1]
        train_qtg(df, train_set) # execute 

    print("--- %s seconds ---" % round((time.time() - start_time),2) ) 

