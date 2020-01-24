#!/usr/bin/python

#!/usr/bin/python

# Version 20200122
# Purpose: Use leave-one-out analysis and the change of AUC-ROC to evaluate importance of each feature. 
# Usage = "QTG_feature_importance.py input_feature_list species_abbreviation"
# feature list: Arabidopsis_features_v5_ortholog.csv for Arabidopsis; rice_features_v5_ortholog.csv for rice ; Setaria_features.csv for Setaria viridis; Sorghum_features.csv for Sorghum 
# species_abbreviation: "AT" for Arabidopsis; "OS" for rice; "SB" for sorghum; "SV" for Setaria
# Usage example: QTG_feature_importance.py Arabidopsis_features_v5_ortholog.csv 'AT'

import numpy as np
import random
import pandas as pd
from sklearn import model_selection
from sklearn.model_selection import StratifiedKFold
from sklearn import svm
from sklearn import metrics
from sklearn import preprocessing
from sklearn import ensemble
from sklearn.metrics import roc_curve, auc, precision_recall_curve, average_precision_score
import sys
from scipy.interpolate import interp1d
import time

start_time = time.time()

def train_qtg_metrics_LOO(df, train_set,out): # out should be a feature name
    if sys.argv[2]=='AT': # species specific opti mparameter
        mol_para=['auto',20]# max_feature, train pos:neg ratio
    if sys.argv[2]=='OS':
        mol_para=['auto',5]
    if sys.argv[2]=='SV':    
        mol_para=['auto',5]
    clf = ensemble.RandomForestClassifier(n_estimators=100,min_samples_split=2,max_features=mol_para[0])
    k_numb_iter=50 ### randomdize k-fold splitting
    all_auc = []

    if  k_numb_iter > 1 : # shuffle on or off based on k_iter
        shuffle_switch=True
    else: 
        shuffle_switch=False
    neg_numb_iter=50 #randomdize the negatives selected

    for n in range(0, k_numb_iter): # these interations are to randomdize k-fold splitting
        skf= StratifiedKFold(n_splits=5,shuffle=shuffle_switch) # the same K split index is used for train and train_L
        roc_auc_mean=[]
        for train_cv, test_cv in skf.split(train_set.drop(['class'], axis=1,inplace=False) ,train_set['class']): #emumerate is to adds a counter
            fold_mean_tpr = np.linspace(1, 0, 101) # set recall
            
            for i in range(0, neg_numb_iter): # iterate n times, randomdize the negatives selected. 
                train_data = train_set.iloc[train_cv] 
                test_data = train_set.iloc[test_cv]
                training_negative = random.sample(list(df[df['class']==0].index), int(len(train_data)*mol_para[1])) # training ratio Os 5 AT 20
                testing_negatives = random.sample(list(df[df['class']==0].index), len(test_cv)*200) # 
                train_data=train_data.append(df.iloc[training_negative]) 
                test_data=test_data.append(df.iloc[testing_negatives])
                train_feature=train_data.drop(['class'], axis=1)
                test_feature=test_data.drop(['class'], axis=1)
                train_feature_L=train_feature.drop(out, axis=1, inplace=False) # drop a feature
                test_feature_L=test_feature.drop(out, axis=1, inplace=False) # drop a feature

                probas_=clf.fit(train_feature, train_data['class']).predict_proba(test_feature) # extract prob of a class (currently 50 test samples)
                probas_L=clf.fit(train_feature_L, train_data['class']).predict_proba(test_feature_L)
                
                fpr, tpr1, thresholds1 = roc_curve(test_data['class'], probas_[:, 1]) # export precision, recall,  associated threshhold at each step
                fpr_L, tpr1_L, thresholds1_L = roc_curve(test_data['class'], probas_L[:, 1]) # export precision, recall,  associated threshhold at each step
                roc_auc_mean.append(auc(fpr_L, tpr1_L)-auc(fpr,tpr1))#calculate AUC reduction
                    
        all_auc.append(np.mean(roc_auc_mean))

    # 
    # average of k iterations
    return ([out,np.mean(all_auc), np.std(all_auc)])
    print ('roc AUC averge %f ; SD is %f' % (np.mean(all_auc), np.std(all_auc))  ) # mean and sd of precision-recall Auc

random.seed(12)
dt = sys.argv[1]
df = pd.read_csv(dt)
df = df.drop(['ID'], axis=1)
df=df.dropna(axis=1,how='all')
df['network_weight'] = preprocessing.scale(df['network_weight'])
train_set = df[df['class']==1]
random.seed(12)
for i in df.columns[:-1]: 
    with open(sys.argv[2]+'_'+'LOO_AUC_reduction.txt','a') as outfile:  #append mode, so delete file before a new run
        a=train_qtg_metrics_LOO(df, train_set,i)
        outfile.write ('{0}\t{1}\t{2}\n'.format(a[0],a[1],a[2]) )

print("--- %s Hrs ---" % round( ((time.time() - start_time)/3600),2) )