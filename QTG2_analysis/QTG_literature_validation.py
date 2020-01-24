
#!/usr/bin/python

# Version 20200122
# Purpose: External validation of models with curated QTL causal genes from literature
# Usage = "QTG_feature_importance.py input_feature_list species_abbreviation"
# feature list: Arabidopsis_features_v5_ortholog.csv for Arabidopsis; rice_features_v5_ortholog.csv for rice ; Setaria_features.csv for Setaria viridis; Sorghum_features.csv for Sorghum 
# species_abbreviation: "AT" for Arabidopsis; "OS" for rice; "SB" for sorghum; "SV" for Setaria
# Usage example: QTG_feature_importance.py Arabidopsis_features_v5_ortholog.csv 'AT'

import numpy as np
import random
import pandas as pd
from sklearn import model_selection
from sklearn import svm
from sklearn import metrics
from sklearn import preprocessing
from sklearn import ensemble
import sys
import time

start_time = time.time()
def train_qtg(df, train_set, Validation_set_i):
    if sys.argv[3]=='AT': # species specific opti mparameter
        mol_para=['auto',20,2]# max_feature, train pos:neg , min sample split
    if sys.argv[3]=='OS':
        mol_para=['auto',5,2] 
    if sys.argv[3]=='SV':
        mol_para=['auto',5,6] # parameter for Setaria
    if sys.argv[3]=='SB':
        mol_para=['auto',5,6] # parameter for sorghum  
    clf = ensemble.RandomForestClassifier(n_estimators=200, min_samples_split=mol_para[2] ,max_features=mol_para[0]) # 50 trees
    neg_inter=5000  # with 5000 repeat, almost no variation

    prediction_list = (len(Validation_set_i))*[0]
    for i in range(0, neg_inter):
                train_data = train_set # positives
                training_negative = random.sample(list(df[df['class']==0].index), len(train_data)*mol_para[1]) # traning neg, pos:neg 5 rice and 20 for At
                train_data=train_data.append(df.iloc[training_negative]) #adding selected negatives to taining dataset
                train_feature=train_data.drop(['class'], axis=1) # drop used columns
                clf.fit(train_feature, train_data['class'])  # model fitting
                validation_feature_all = Validation_set_i.drop(['class'], axis=1) 
                validation_pred = clf.predict_proba(validation_feature_all)[:,1] 
                prediction_list+= validation_pred
    prediction_list_all=np.array(prediction_list)#required format for next step 
    Ranks=prediction_list_all.argsort()[::-1].argsort()
    causal_index=Validation_set_i.index[Validation_set_i['class']==1].tolist() # extract causal index
    causal_rank=Ranks[causal_index]
    prediction_list_freq= [l/(neg_inter) for l in prediction_list]
    with open('identification_frequency.csv','w') as indenti_f: 
        for i in range(len(prediction_list_freq)):
            indenti_f.write (Validation_set_ID_uni[i]+','+str(prediction_list_freq[i])+'\n')

    for i in causal_index:
        print('Causal gene {0} rank {1} ({2}%), with {3} interation'.format(Validation_set_ID_uni[i],Ranks[i]+1,int(((Ranks[i])/original_length)*100),neg_inter)) # ranking of each validation sample

random.seed(12)
dt = sys.argv[1]
print(dt)
df = pd.read_csv(dt)
df=df.dropna(axis=1,how='all')



dt1= sys.argv[2]
Validation_set_ID=pd.read_csv(dt1,names=['ID','class'],header=None)

ID_list=Validation_set_ID[Validation_set_ID['class']==1]['ID'].tolist()# this to extract causal genes in validation set and remove the genes from training set 
for i in ID_list: 
    index_interm=df[df['ID']==1].index
    df.loc[index_interm,'class']=0

original_length= len(Validation_set_ID) # original number of gene on QTL  

Validation_set=pd.DataFrame()   # not all genes will be in feature list, so the gene have to be removed. 
for i in range(len(Validation_set_ID)):   
    Validation_set=Validation_set.append (df[df.ID==Validation_set_ID.ID[i]])

for i in range(len(Validation_set)):
    Validation_set.iloc[i,len(Validation_set.columns)-1]= int(Validation_set_ID[Validation_set_ID.ID==Validation_set.iloc[i,0]]['class'])

Validation_set_ID_uni=list(Validation_set.ID) 
trimed_length=len(set(Validation_set_ID_uni)) # length of genes if only gene annottions shared by RAP and MSU are considered
dropped_length=original_length-trimed_length


print('Number of genes: '+str(original_length)  ) 

df = df.drop(['ID'], axis=1) # remove the ID column
Validation_set= Validation_set.drop(['ID'], axis=1) # remove the ID column
Validation_set_i=Validation_set.reset_index(drop=True)

## assign validation_set and train_set

train_set = df[df['class']==1]


train_qtg(df, train_set, Validation_set_i)

print("--- %s seconds ---" % round((time.time() - start_time),2) )