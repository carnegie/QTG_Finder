#!/usr/bin/python

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
    if sys.argv[2]=='AT': # species specific opti mparameter
        mol_para=['auto',20]# max_feature, train pos:neg ratio
    if sys.argv[2]=='OS':
        mol_para=['auto',5]
    clf = ensemble.RandomForestClassifier(n_estimators=200,min_samples_split=2,max_features=mol_para[0])
    # fout = open("P_R.txt", 'w')
    froc =open(sys.argv[2]+"_"+"ROC.txt", 'w')
    # fimp = open('feature_importance.txt', 'w')
    k_numb_iter=50 ### randomdize k-fold splitting
    all_mean_prec = [float()]*k_numb_iter
    all_mean_tpr1 = [float()]*k_numb_iter
    P_R_all_auc = []
    ROC_all_auc=[]
    importance=[]
    if  k_numb_iter > 1 : # re-shuffle on or off based on k_iter
        shuffle_switch=True
    else: 
        shuffle_switch=False
    neg_numb_iter=100 # iteration for randomly select negatives 

    for n in range(0, k_numb_iter): # these interations are to randomdize k-fold splitting

        skf=KFold(len(train_set), n_folds=5,shuffle=shuffle_switch)
        P_R_auc_mean=[]
        ROC_auc_mean=[]
        fold_mean_prec = [0.0]*101 
        fold_mean_tpr1 = [0.0]*101 
        for m, (train_cv, test_cv) in enumerate(skf): #emumerate is to adds a counter
    
            fold_mean_tpr = np.linspace(0, 1, 101) # set recall
            fold_mean_fpr = np.linspace(0, 1, 101) # set fpr grid for ROC

            for i in range(0, neg_numb_iter): # iterate n times, randomdize the negatives selected. 
                train_data = train_set.iloc[train_cv] 
                test_data = train_set.iloc[test_cv] 
                training_negative = random.sample(list(df[df['class']==0].index), int(len(train_data)*mol_para[1])) # train neg:pos=1:1
                testing_negatives = random.sample(list(df[df['class']==0].index), len(test_cv)*200) # test negatives are randomly selected, neg:pos=200:1
                train_data=train_data.append(df.iloc[training_negative]) 
                test_data=test_data.append(df.iloc[testing_negatives])
                train_feature=train_data.drop(['class'], axis=1)
                test_feature=test_data.drop(['class'], axis=1)
                probas_=clf.fit(train_feature, train_data['class']).predict_proba(test_feature) # extract prob of a class (currently 50 test samples)
                importance.append(clf.feature_importances_) # extract importance of features (currently 28). 
                prec, tpr, thresholds = precision_recall_curve(test_data['class'], probas_[:, 1]) # export precision, recall,  associated threshhold at each step
                fpr, tpr1, thresholds1 = roc_curve(test_data['class'], probas_[:, 1])
                ROC_auc_mean.append(auc(fpr, tpr1))
                P_R_auc_mean.append(auc(tpr, prec)) #calculate AUC at here is more accurate
                fun=interp1d( tpr,  prec) # use interpolation to estimate the PR curve
                fold_mean_prec += fun(fold_mean_tpr) # cal standardized precision with tpr grid
                fun2=interp1d(fpr, tpr1)
                fold_mean_tpr1 += fun2(fold_mean_fpr) # cal standardized fpr with tpr1 grid

        fold_mean_prec /= 5*neg_numb_iter #  the average of 5 cr runs, neg iterations  
        fold_mean_tpr1 /= 5*neg_numb_iter
        all_mean_prec[n] = fold_mean_prec # accumulative cr results # need to make this a dataframe to calc SD
        all_mean_tpr1[n]= fold_mean_tpr1
        P_R_all_auc.append(np.mean(P_R_auc_mean))
        ROC_all_auc.append(np.mean(ROC_auc_mean))

    # sd of iterations
    all_mean_prec=np.array(all_mean_prec)
    average_prec = [float()]*101 # 101 is grid number,see fold_mean_tpr
    sd_prec = [float()]*101
    for i in range(101):
        average_prec[i]=np.mean(all_mean_prec[:,i])
        sd_prec[i]=np.std(all_mean_prec[:,i])

    
    all_mean_tpr1=np.array(all_mean_tpr1)
    average_tpr1 = [float()]*101 # 101 is grid number,see fold_mean_tpr
    sd_tpr1 = [float()]*101
    for i in range(101):
        average_tpr1[i]=np.mean(all_mean_tpr1[:,i])
        sd_tpr1[i]=np.std(all_mean_tpr1[:,i])


    # average of k iterations
    
    # print ('P-R AUC averge %f ; SD is %f' % (np.mean(P_R_all_auc), np.std(P_R_all_auc))  ) # mean and sd of precision-recall Auc
    print ('ROC AUC averge %f ; SD is %f' % (np.mean(ROC_all_auc), np.std(ROC_all_auc))  ) # mean and sd of roc Auc

    mean_importance = [0.0] * len(importance[0]) # creat a empty vector for mean importance 
    std_importance = [0.0] * len(importance[0])
    for i in range(len(importance[0])): # number of features
        imp = []
        for j in range(len(importance)): # number of repeat to be averaged
            mean_importance[i] += importance[j][i]
            imp.append(importance[j][i])
        mean_importance[i] = np.mean(imp)  
        std_importance[i] = np.std(imp)

    indices = np.argsort(mean_importance)[::-1] # return index for sorting, descending
    #for i in range(len(mean_importance)):
    #    fimp.write("%d. feature %s (%f, %f)\n" % (i + 1, df.columns[indices[i]], mean_importance[indices[i]], std_importance[indices[i]]))

    #for i in range(0, len(fold_mean_tpr)):
    #    fout.write('{0}\t{1}\t{2}\n'.format(fold_mean_tpr[i], average_prec[i], sd_prec[i]))
    
    for i in range(0, len(fold_mean_fpr)):
        froc.write('{0}\t{1}\t{2}\n'.format(fold_mean_fpr[i], average_tpr1[i], sd_tpr1[i]))


dt = sys.argv[1]
df = pd.read_csv(dt)
df = df.drop(['ID'], axis=1)
df=df.dropna(axis=1,how='all')


df['network_weight'] = preprocessing.scale(df['network_weight'])

train_set = df[df['class']==1]

train_qtg_metrics(df, train_set)

print("--- %s seconds ---" % round((time.time() - start_time),2) )