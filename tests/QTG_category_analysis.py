#!/usr/bin/python

# Version 20181106
# Purpose: Use leave-one-out analysis to evaluate if the algrithm prefer certain categories of causal genes over other categories
# Usage = "QTG_literature_validation.py input_feature_list species_abbreviation"
# species_abbreviation: "AT" for Arabidopsis; "OS" for rice
# Usage example: QTG_literature_validation.py Arabidopsis_features-v3.05.txt 'AT'  


import numpy as np
import random
import pandas as pd
from sklearn import cross_validation
from sklearn.cross_validation import KFold
from sklearn import svm
from sklearn import metrics
from sklearn import preprocessing
from sklearn import ensemble

import sys

def train_qtg(df, train_set, validation_set):
    if sys.argv[3]=='AT': 
        mol_para=['auto',20]# optimum parameters for Arabidopsis
    if sys.argv[3]=='OS':
        mol_para=['auto',5] # optimum parameters for rice
    clf = ensemble.RandomForestClassifier(n_estimators=200, min_samples_split=2,max_features=mol_para[0]) # 200 trees
    fout = open(sys.argv[3]+"category_rank.txt", 'w') # output file
    prediction_list = [[] for x in range(len(validation_set))]
    neg_inter=5000 #

    for inx in range(len(validation_set)):
        for i in range(0, neg_inter):
            train_data =train_set.drop(validation_set.index[inx] , inplace=False) # remove one causal gene from the training set 
            training_negative = random.sample(list(df[df['class']==0].index), len(train_data)*mol_para[1]) # # randomly select negatives from genome genes 
            train_data=train_data.append(df.iloc[training_negative]) 
            train_feature=train_data.drop(['class'], axis=1) 
            clf.fit(train_feature, train_data['class'])  # model fitting
            m = validation_set.index[inx] # pull out index of the causal for validation
            upperbound=m+100 # 100 genes after the causal gene 
            lowerbound=m-100 # 100 genes before the causal gene 
            validation_set_all=df.iloc[lowerbound:upperbound] # a list with a causal gene and 200 genes nearby the causal gene 
            validation_feature_all = validation_set_all.drop(['class'], axis=1) 
            validation_pred =clf.predict_proba (validation_feature_all)[:,1] # apply model to the causal gene and 200 genes nearby
            if i == 0:
                prediction_list[inx] = validation_pred/neg_inter 
            else:
                prediction_list[inx] = [sum(x) for x in zip(prediction_list[inx], validation_pred/neg_inter)]# combine prediction results 
    order_cv = [0]*len(validation_set)
    for inx in range(len(validation_set)):
        order=sorted(range(len(prediction_list[inx])), key = lambda k:prediction_list[inx][k], reverse=True) # rank all genes based on times of being identified as a causal gene
        order_cv[inx] += order.index(100) # assign the ranks,index 100 is the position casual gene  

    print(order_cv) # ranking of each validation sample
    print(len(order_cv))
    print(len(validation_set_copy))
    print(validation_set_copy.ID)

    for i in range(len(validation_set_copy)):
        fout.write (validation_set_copy.ID[i]+'\t'+str(order_cv[i])+'\n') # out putfile with the rank of each causal gene among 200 nearby genes
    fout.write('\n'.join([str(i) for i in order_cv]))
    fout.close()
    
dt = sys.argv[1]
df = pd.read_csv(dt)

dt1= sys.argv[2]
Validation_set_ID=pd.read_csv(dt1)

Vali_index=[]
for i in range(len(Validation_set_ID)):
    Vali_index.extend(df[df.ID==Validation_set_ID.ID[i]].index.tolist())
    
validation_set =df.iloc[Vali_index]
validation_set_copy= validation_set.copy().reset_index()

validation_set=validation_set.drop(['ID'], axis=1)
validation_set['network_weight'] = preprocessing.scale(validation_set['network_weight'])
validation_set=validation_set.dropna(axis=1)

df = df.drop(['ID'], axis=1) 
df['network_weight'] = preprocessing.scale(df['network_weight'])
df = df.dropna(axis=1)
train_set = df[df['class']==1]

train_qtg(df, train_set, validation_set)