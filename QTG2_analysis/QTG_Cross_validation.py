#!/usr/bin/python

# Version 20200113
# Purpose: Use cross validaiton to calculate the AUC-ROC 
# Usage = "QTG_Cross_validation.py input_feature_list species_abbreviation"
# feature list: Arabidopsis_features_v5_ortholog.csv for Arabidopsis; rice_features_v5_ortholog.csv for rice ; Setaria_features.csv for Setaria viridis; Sorghum_features.csv for Sorghum 
# species_abbreviation: "AT" for Arabidopsis; "OS" for rice; "SB" for sorghum; "SV" for Setaria
# Usage example: QTG_Cross_validation.py Arabidopsis_features_v5_ortholog.csv 'AT'

import numpy as np
import random
import pandas as pd
from sklearn import model_selection 
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import confusion_matrix
from sklearn.metrics import f1_score
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
        mol_para=['auto',20,2]# parameter for Arabidopsis
    if sys.argv[2]=='OS':
        mol_para=['auto',5,2] # parameter for rice
    if sys.argv[2]=='SV':
        mol_para=['auto',5,6] # parameter for Setaria
    if sys.argv[2]=='SB':
        mol_para=['auto',5,6] # parameter for Sorghum
    clf = ensemble.RandomForestClassifier(n_estimators=100,min_samples_split=mol_para[2],max_features=mol_para[0]) # random forest parameters
    froc =open(sys.argv[2]+"_"+"ROC.txt", 'w') # output file for ploting ROC curve
    fPR =open(sys.argv[2]+"_"+"PR.txt", 'w') # output file for ploting ROC curve
    k_numb_iter=50 #  interations for randomly re-spliting causal gene list, 30
    all_mean_prec = [float()]*k_numb_iter
    all_mean_tpr1 = [float()]*k_numb_iter
    P_R_all_auc = []
    ROC_all_auc=[]
    importance=[]
    F1_all=[]
    if  k_numb_iter > 1 : # re-shuffle on or off based on k_iter
        shuffle_switch=True
    else: 
        shuffle_switch=False
    neg_numb_iter=30 # iteration for randomly select negatives and re-train model, 50
    for n in range(0, k_numb_iter): #  re-spliting causal gene list for cross-validation
        #skf=KFold(len(train_set), n_folds=5,shuffle=shuffle_switch) # five fold cross-validatoin
        skf= StratifiedKFold(n_splits=5,shuffle=shuffle_switch)
        P_R_auc_mean=[]
        ROC_auc_mean=[]
        confusion_matrix_sum=np.array([0,0,0,0])
        F1_sum=0
        fold_mean_prec = [0.0]*101 
        fold_mean_tpr1 = [0.0]*101 
        temp_actual_label=[]
        temp_predicted_class=[]
        for train_cv, test_cv in skf.split(train_set.drop(['class'], axis=1,inplace=False) ,train_set['class']): 
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
                predicted_class=clf.fit(train_feature, train_data['class']).predict(test_feature) # for confusion matrix and F1
                #importance.append(clf.feature_importances_) # feature importance (Gini importance), not used 
                prec, tpr, thresholds = precision_recall_curve(test_data['class'], probas_[:, 1]) # precision, recall,  associated threshhold 
                fpr, tpr1, thresholds1 = roc_curve(test_data['class'], probas_[:, 1]) # false positve rate, true positiverate, and treshhold
                tn, fp, fn, tp = confusion_matrix(test_data['class'], predicted_class).ravel()
                temp_actual_label.extend(test_data['class'])
                temp_predicted_class.extend( predicted_class)
                confusion_matrix_sum+= np.array([tn, fp, fn, tp])
                ROC_auc_mean.append(auc(fpr, tpr1)) # AUC-ROC 
                P_R_auc_mean.append(auc(tpr, prec)) # AUC-precision_recall
                fun2=interp1d(fpr, tpr1) 
                fold_mean_tpr1 += fun2(fold_mean_fpr) # fit fpr with tpr1 grid
                fun3=interp1d (tpr, prec)
                fold_mean_prec+=fun3(fold_mean_fpr) 
        
        F1_score=f1_score(temp_actual_label, temp_predicted_class)
        F1_all.append(F1_score)
        fold_mean_tpr1 /= 5*neg_numb_iter #  the average of cross validation iterations  
        all_mean_tpr1[n]= fold_mean_tpr1 
        ROC_all_auc.append(np.mean(ROC_auc_mean))
        P_R_all_auc.append(np.mean(P_R_auc_mean))
        fold_mean_prec/=5*neg_numb_iter
        all_mean_prec[n]= fold_mean_prec

    all_mean_tpr1=np.array(all_mean_tpr1)
    average_tpr1 = [float()]*101 # 101 is grid number,see fold_mean_tpr
    sd_tpr1 = [float()]*101
    for i in range(101):
        average_tpr1[i]=np.mean(all_mean_tpr1[:,i]) # average of true positive rate at each grid point
        sd_tpr1[i]=np.std(all_mean_tpr1[:,i])
    
    print ('ROC AUC averge %f ; SD is %f' % (np.mean(ROC_all_auc), np.std(ROC_all_auc))  ) # mean and standard deviation of AUC-ROC    
    print ('precision recall AUC averge %f ; SD is %f' % (np.mean(P_R_all_auc), np.std(P_R_all_auc))  )
    for i in range(0, len(fold_mean_fpr)):
        froc.write('{0}\t{1}\t{2}\n'.format(fold_mean_fpr[i], average_tpr1[i], sd_tpr1[i]))
    
    all_mean_prec=np.array(all_mean_prec)
    average_prec = [float()]*101
    sd_prec = [float()]*101
    for i in range(101):
        average_prec[i]=np.mean(all_mean_prec[:,i])
        sd_prec[i]=np.std(all_mean_prec[:,i])    
    
    for i in range(0, len(fold_mean_fpr)):
        fPR.write ('{0}\t{1}\t{2}\n'.format(fold_mean_fpr[i], average_prec[i], sd_prec[i]))

    confusion_matrix_avg=confusion_matrix_sum/(5*neg_numb_iter)
    F1_avg=F1_sum/(5*neg_numb_iter)
    print('\nConfusion matrix (TN, FP, FN, TP, average of all iterations and K-folds):\n',confusion_matrix_avg)
    print('\nF1 score average:\n',np.mean(F1_all),'F1 score sd:\n',np.std(F1_all))

random.seed(12)
dt = sys.argv[1] # input feature list
df = pd.read_csv(dt)
df = df.drop(['ID'], axis=1)
df=df.dropna(axis=1,how='all')

#df['network_weight'] = preprocessing.scale(df['network_weight'])
train_set = df[df['class']==1]
train_qtg_metrics(df, train_set)

print("--- %s seconds ---" % round((time.time() - start_time),2) )