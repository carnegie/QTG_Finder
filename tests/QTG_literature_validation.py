#!/usr/bin/python

# Version 20181106
# Purpose: Use leave-one-out analysis and the change of AUC-ROC to evaluate importance of each feature. Note: This script is for testing causal genes in a single QTL. Batch files 'AT_lit_validation.sh' or 'OS_lit_validation.sh' are abaliable to analyze multiple QTLs.  
# Usage = "QTG_literature_validation.py feature_list species_abbreviation QTL_gene_list_flag"
# feature list: use Arabidopsis_features_v4.csv for Arabidopsis; use rice_features_v2.csv for rice 
# species_abbreviation: "AT" for Arabidopsis; "OS" for rice
# QTL_gene_list_flag: a list of genes in the QTL. Causal genes are labeled with 1 in the second column. Other genome genes are labeled as 0. examples can be found in 'tests/input/AT' or 'tests/input/OS/'
# Usage example: QTG_literature_validation.py Arabidopsis_features_v4.csv 'AT' ./input/AT/Huang_2012.csv  

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
import time

start_time = time.time()
def train_qtg(df, train_set, Validation_set_i):
    if sys.argv[3]=='AT': 
        mol_para=[9,20]# parameters for Arabidopsis
    if sys.argv[3]=='OS':
        mol_para=[9,5] # parameters for rice  
    clf = ensemble.RandomForestClassifier(n_estimators=100, min_samples_split=2,max_features=mol_para[0]) 
    neg_inter=5000  # interrations for randomly selecting negatives from genome and re-training models 

    prediction_list = (len(Validation_set_i))*[0]
    for i in range(0, neg_inter):
                train_data = train_set 
                training_negative = random.sample(list(df[df['class']==0].index), len(train_data)*mol_para[1]) # randomly select negatives from genome genes 
                train_data=train_data.append(df.iloc[training_negative]) 
                train_feature=train_data.drop(['class'], axis=1) 
                clf.fit(train_feature, train_data['class'])  # model fitting
                validation_feature_all = Validation_set_i.drop(['class'], axis=1) 
                validation_pred = clf.predict_proba(validation_feature_all)[:,1] 
                prediction_list+= validation_pred
    prediction_list_all=np.array(prediction_list)
    Ranks=prediction_list_all.argsort()[::-1].argsort()
    causal_index=Validation_set_i.index[Validation_set_i['class']==1].tolist() # extract causal index
    causal_rank=Ranks[causal_index]
    prediction_list_freq= [l/(neg_inter) for l in prediction_list]
    #with open('identification_frequency.csv','w') as indenti_f: # frequency of beging identified as a causal gene, used for trouble shooting  
    #    for i in range(len(prediction_list_freq)):
    #        indenti_f.write (Validation_set_ID_uni[i]+','+str(prediction_list_freq[i])+'\n')

    for i in causal_index:
        print('Causal gene {0} rank {1} ({2}%), with {3} interation'.format(Validation_set_ID_uni[i],Ranks[i]+1,int(((Ranks[i])/original_length)*100),neg_inter)) # ranking of each causal gene in QTL

dt = sys.argv[1]
df = pd.read_csv(dt)
df=df.dropna(axis=1,how='all')

df['network_weight'] = preprocessing.scale(df['network_weight'])

dt1= sys.argv[2]
Validation_set_ID=pd.read_csv(dt1,names=['ID','class'],header=None)

original_length= len(Validation_set_ID) # number of gene on the QTL  

Validation_set=pd.DataFrame()   # exclude genes not in the feature list
for i in range(len(Validation_set_ID)):   
    Validation_set=Validation_set.append (df[df.ID==Validation_set_ID.ID[i]])
for i in range(len(Validation_set)):
    Validation_set.iloc[i,len(Validation_set.columns)-1]= int(Validation_set_ID[Validation_set_ID.ID==Validation_set.iloc[i,0]]['class'])

Validation_set_ID_uni=list(Validation_set.ID) 

print('Number of genes: '+str(original_length)  ) 

## assign validation_set and train_set
df = df.drop(['ID'], axis=1) 
Validation_set=  Validation_set.drop(['ID'], axis=1) 
Validation_set_i=Validation_set.reset_index(drop=True)

train_set = df[df['class']==1]
train_qtg(df, train_set, Validation_set_i)

print("--- %s seconds ---" % round((time.time() - start_time),2) )