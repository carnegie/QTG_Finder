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
import sys
import time
import itertools as it

start_time = time.time()

def train_qtg(df, train_set, Validation_set_i):
    clf = ensemble.RandomForestClassifier(n_estimators=200, min_samples_split=2) # Random forest parameters
    neg_inter=5000  # interrations of randomly select negatives from genome and re-train models 
    prediction_list = (len(Validation_set_i))*[0]
    for i in range(0, neg_inter):
                train_data = train_set # positives used for training 
                training_negative = random.sample(list(df[df['class']==0].index), len(train_data)*20) # randomly select negatives from genome genes
                train_data=train_data.append(df.iloc[training_negative]) 
                train_feature=train_data.drop(['class'], axis=1) 
                clf.fit(train_feature, train_data['class'])  # model fitting
                validation_feature_all = Validation_set_i.drop(['class'], axis=1) 
                validation_pred = clf.predict_proba(validation_feature_all)[:,1] # prediction based on models
                prediction_list+= validation_pred # combine prediction results
    prediction_list_all=np.array(prediction_list) 
    Ranks=prediction_list_all.argsort()[::-1].argsort() # rank gene based on times of being predicted as a causal gene
    #causal_index=Validation_set_i.index[Validation_set_i['class']==1].tolist() # extract causal index
    #causal_rank=Ranks[causal_index]
    prediction_list_freq= [l/(neg_inter) for l in prediction_list] # cal freq identified as causal 
    rank_list_df= pd.DataFrame( {'ID':Validation_set_ID_uni, 'freq': prediction_list_freq })
    
    for i in gene_ex:  # not included this step if want to use known causal as pos control
        rank_list_df= rank_list_df.append({'ID': i, 'freq':0 }, ignore_index=True) # put excluded gene back with a freq of 0, need to manually remove known causal gene from list if necessary  

    rank_list_df['Rank'] = rank_list_df['freq'].rank(ascending=0,method='average') # creat rank column
    rank_list_df_sorted=rank_list_df.sort_values(by=['Rank'],ascending =True) #sort df
    rank_list_df_sorted=rank_list_df_sorted.reset_index(drop= True)

    with open('QTL_gene_rank.csv','a') as indenti_f: 
        indenti_f.write ('//'+'\n'+QTL_name+'\n')
        indenti_f.write ('ID'+','+'Rank_in_a_QTL'+','+'Frequency'+'\n')
        for i in range(len(rank_list_df_sorted)):            
            indenti_f.write (rank_list_df_sorted['ID'][i]+','+str(rank_list_df_sorted['Rank'][i])+','+ str(rank_list_df_sorted['freq'][i]) +'\n')
            #indenti_f.write (rank_list_df['ID'][i]+','+str(rank_list_df['Rank'][i])+','+ str(rank_list_df['freq'][i])+'\n')  # alternative, including freq score

dt = sys.argv[1] # used to extract causal gene list

#### check know QTL gene on QTLs, remove, reducant ID run

#dt1= sys.argv[2]
#Validation_set_ID=pd.read_csv(dt1,names=['ID','class'],header=None)

with open(sys.argv[2],'r') as f:  
        for key,group in it.groupby(f,lambda line: line.startswith('//') ):
            if not key:
                df = pd.read_csv(dt)
                df=df.dropna(axis=1,how='all')
                df['network_weight'] = preprocessing.scale(df['network_weight'])
                group = list(group)
                group=[i.strip('\n') for i in group]
                QTL_name=group[0]
                print('QTL: '+QTL_name)
                genes_in_QTL=group[1:]
                original_length= len(genes_in_QTL) # original number of gene on QTL  
                Validation_set=pd.DataFrame()   # not all genes will be in feature list, so the gene have to be removed. 
                for i in range(len(genes_in_QTL)):   
                    Validation_set=Validation_set.append (df[df.ID==genes_in_QTL[i]])
                df = df.drop(['ID'], axis=1) # remove the ID column # repeat dropping is a problem?
                ind_for_exclusion=[]
                for t in range(len(Validation_set.index)): 
                    if Validation_set['class'][Validation_set.index[t]]==1:
                        ind_for_exclusion.append(t)
                        print ('Known QTL causal gene excluded: ')
                        print(Validation_set['ID'][Validation_set.index[t]])  # known causal in QTL
                #print(Validation_set.index[ind_for_exclusion]) # index for exclusion
                #print(len(Validation_set))
                Validation_set=Validation_set.drop(Validation_set.index[ind_for_exclusion])  # take out known causal from data,
                #print(len(Validation_set))
                Validation_set_ID_uni=list(Validation_set.ID) # ID for validation set, order not change
                trimed_length=len(set(Validation_set_ID_uni)) # length of genes if only gene annottions shared by RAP and MSU are considered
                print('Number of genes: '+str(original_length)  )  #sub: if only consider gene annottions shared by RAP and MSU, also excluding known causal
                gene_ex=(set(genes_in_QTL)-set(Validation_set_ID_uni)) # gene not shared in  RAp and MSU and known causal gene 
                Validation_set=  Validation_set.drop(['ID'], axis=1) # remove the ID column
                Validation_set_i=Validation_set.reset_index(drop=True)
                train_set = df[df['class']==1]## assign  train_set
                #print(len(Validation_set_i))
                train_qtg(df, train_set, Validation_set_i) # run

print("--- %s seconds ---" % round((time.time() - start_time),2) )

