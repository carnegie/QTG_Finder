#!/usr/bin/python

# Version 20181106
# Purpose: Use leave-one-out analysis and the change of AUC-ROC to evaluate importance of each feature. 
# Usage = "QTG_feature_importance.py input_feature_list species_abbreviation"
# feature list: use Arabidopsis_features-v3.05_n.csv for Arabidopsis; use rice_features_v1.3.11_n.csv for rice 
# species_abbreviation: "AT" for Arabidopsis; "OS" for rice
# Usage example: QTG_feature_importance.py Arabidopsis_features-v3.05_n.csv 'AT'

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

def train_qtg_metrics_LOO(df, train_set,out): 
    if sys.argv[2]=='AT': 
        mol_para=['auto',20]# parameters for Arabidopsis
    if sys.argv[2]=='OS':
        mol_para=['auto',5]# parameters for Arabidopsis
    clf = ensemble.RandomForestClassifier(n_estimators=200,min_samples_split=2,max_features=mol_para[0]) # Random forest parameters
    k_numb_iter=50 #  interations for randomly re-spliting causal gene list
    all_auc = []

    if  k_numb_iter > 1 : # re-shuffle on or off based on k_iter
        shuffle_switch=True
    else: 
        shuffle_switch=False

    neg_numb_iter=100 #interrations for randomly selecting negatives from genome and re-training models 
    for n in range(0, k_numb_iter): 
        skf=KFold(len(train_set), n_folds=5,shuffle=shuffle_switch) # five fold cross-validatoin
        roc_auc_mean=[]
        for m, (train_cv, test_cv) in enumerate(skf): 
            fold_mean_tpr = np.linspace(1, 0, 101) 
            for i in range(0, neg_numb_iter): # terate for randomly selecting negatives. 
                train_data = train_set.iloc[train_cv] 
                test_data = train_set.iloc[test_cv]
                training_negative = random.sample(list(df[df['class']==0].index), int(len(train_data)*mol_para[1])) 
                testing_negatives = random.sample(list(df[df['class']==0].index), len(test_cv)*200) # 
                train_data=train_data.append(df.iloc[training_negative]) 
                test_data=test_data.append(df.iloc[testing_negatives])
                train_feature=train_data.drop(['class'], axis=1)
                test_feature=test_data.drop(['class'], axis=1)
                train_feature_L=train_feature.drop(out, axis=1, inplace=False) # leave one feature out
                test_feature_L=test_feature.drop(out, axis=1, inplace=False) # leave one feature out

                probas_=clf.fit(train_feature, train_data['class']).predict_proba(test_feature) # the predictions of full model 
                probas_L=clf.fit(train_feature_L, train_data['class']).predict_proba(test_feature_L) # the predictions of leave-one-out model 
                
                fpr, tpr1, thresholds1 = roc_curve(test_data['class'], probas_[:, 1]) # export precision, recall,  associated threshhold of the full model 
                fpr_L, tpr1_L, thresholds1_L = roc_curve(test_data['class'], probas_L[:, 1]) # export precision, recall,  associated threshhold of the leave-one-out model
                roc_auc_mean.append(auc(fpr_L, tpr1_L)-auc(fpr,tpr1)) #calculate the change of AUC-ROC (delta AUC-ROC)
                    
        all_auc.append(np.mean(roc_auc_mean)) 

    return ([out,np.mean(all_auc), np.std(all_auc)])# average and standarad deviation of delta AUC-ROC for each feature
    print ('AUC-ROC averge %f ; SD is %f' % (np.mean(all_auc), np.std(all_auc))  ) 

dt = sys.argv[1]
df = pd.read_csv(dt)
df = df.drop(['ID'], axis=1)
df=df.dropna(axis=1,how='all')
df['network_weight'] = preprocessing.scale(df['network_weight'])
train_set = df[df['class']==1]

for i in df.columns[:-1]: 
    with open(sys.argv[2]+'_'+'LOO_AUC_reduction.txt','a') as outfile:  #appending mode, so delete file before a new run
        a=train_qtg_metrics_LOO(df, train_set,i)
        outfile.write ('{0}\t{1}\t{2}\n'.format(a[0],a[1],a[2]) ) # output feature importance results

print("--- %s Hrs ---" % round( ((time.time() - start_time)/3600),2) )