#!/usr/bin/python

# Version 20190307
# Purpose: The AUC-ROC from cross validaiton was used to optimize model parameters
# Usage = "QTG_Cross_validation.py input_feature_list species_abbreviation"
# feature list: use Arabidopsis_features_v4 for Arabidopsis; use rice_features_v2.csv for rice 
# species_abbreviation: "AT" for Arabidopsis; "OS" for rice
# Usage example: QTG_Cross_validation.py Arabidopsis_features_v4 'AT'

import numpy as np
import random
import pandas as pd
from sklearn import cross_validation
from sklearn.cross_validation import KFold
from sklearn import svm
from sklearn import metrics
from sklearn import preprocessing
from sklearn import ensemble
from sklearn.metrics import roc_curve, auc, precision_recall_curve, average_precision_score
import sys
from scipy.interpolate import interp1d
import time

start_time = time.time()

def train_qtg_metrics(df,train_set):
    if sys.argv[2]=='AT': 
        mol_para=[9,20]# parameter for Arabidopsis
    if sys.argv[2]=='OS':
        mol_para=[9,5]# parameter for rice
    k_numb_iter=50 #  interations for randomly re-spliting causal gene list
    clf = ensemble.RandomForestClassifier(n_estimators=100,min_samples_split=2,max_features=mol_para[0]) # random forest parameters
    froc =open(sys.argv[2]+"_"+"ROC.txt", 'w') # output file for ploting ROC curve
    all_mean_prec = [float()]*k_numb_iter
    all_mean_tpr1 = [float()]*k_numb_iter
    P_R_all_auc = []
    ROC_all_auc=[]
    importance=[]
    if  k_numb_iter > 1 : # re-shuffle on or off based on k_iter
        shuffle_switch=True
    else: 
        shuffle_switch=False
    neg_numb_iter=50 # iteration for randomly select negatives and re-train model
    for n in range(0, k_numb_iter): #  re-spliting causal gene list for cross-validation
        skf=KFold(len(train_set), n_folds=5,shuffle=shuffle_switch) # five fold cross-validatoin
        P_R_auc_mean=[]
        ROC_auc_mean=[]
        fold_mean_prec = [0.0]*101 
        fold_mean_tpr1 = [0.0]*101 
        for m, (train_cv, test_cv) in enumerate(skf): 
            fold_mean_tpr = np.linspace(0, 1, 101) # set recall
            fold_mean_fpr = np.linspace(0, 1, 101) # set fpr grid for ROC
            for i in range(0, neg_numb_iter): # iterate for randomly selecting negatives. 
                train_data = train_set.iloc[train_cv] 
                test_data = train_set.iloc[test_cv] 
                training_negative = random.sample(list(df[df['class']==0].index), int(len(train_data)*mol_para[1])) # randomly select training negatives from genome genes
                testing_negatives = random.sample(list(df[df['class']==0].index), len(test_cv)*200) # randomly select testing negatives from genome genes, negative:positive=200:1
                train_data=train_data.append(df.iloc[training_negative]) 
                test_data=test_data.append(df.iloc[testing_negatives])
                train_feature=train_data.drop(['class'], axis=1)
                test_feature=test_data.drop(['class'], axis=1)
                probas_=clf.fit(train_feature, train_data['class']).predict_proba(test_feature) #  Probability of falling into a class 
                importance.append(clf.feature_importances_) # feature importance (Gini importance), not used 
                prec, tpr, thresholds = precision_recall_curve(test_data['class'], probas_[:, 1]) # precision, recall,  associated threshhold 
                fpr, tpr1, thresholds1 = roc_curve(test_data['class'], probas_[:, 1]) # false positve rate, true positiverate, and treshhold
                ROC_auc_mean.append(auc(fpr, tpr1)) # AUC-ROC 
                P_R_auc_mean.append(auc(tpr, prec)) # AUC-precision_recall
                fun2=interp1d(fpr, tpr1)
                fold_mean_tpr1 += fun2(fold_mean_fpr) # fit fpr with tpr1 grid

        fold_mean_tpr1 /= 5*neg_numb_iter #  the average of cross validation iterations  
        all_mean_tpr1[n]= fold_mean_tpr1 
        ROC_all_auc.append(np.mean(ROC_auc_mean))

    all_mean_tpr1=np.array(all_mean_tpr1)
    average_tpr1 = [float()]*101 # 101 is grid number,see fold_mean_tpr
    sd_tpr1 = [float()]*101
    for i in range(101):
        average_tpr1[i]=np.mean(all_mean_tpr1[:,i]) # average of true positive rate at each grid point
        sd_tpr1[i]=np.std(all_mean_tpr1[:,i])
    
    print ('ROC AUC averge %f ; SD is %f' % (np.mean(ROC_all_auc), np.std(ROC_all_auc))  ) # mean and standard deviation of AUC-ROC    
    for i in range(0, len(fold_mean_fpr)):
        froc.write('{0}\t{1}\t{2}\n'.format(fold_mean_fpr[i], average_tpr1[i], sd_tpr1[i]))

dt = sys.argv[1] # input feature list
df = pd.read_csv(dt)
df = df.drop(['ID'], axis=1)
df=df.dropna(axis=1,how='all')

df['network_weight'] = preprocessing.scale(df['network_weight'])
train_set = df[df['class']==1]
train_qtg_metrics(df, train_set)

print("--- %s seconds ---" % round((time.time() - start_time),2) )