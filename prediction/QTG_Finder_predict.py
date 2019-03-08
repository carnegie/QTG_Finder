#!/usr/bin/python

# Version 1.1
# last update: 20190307
# Purpose: Use for ranking causal genes in QTL regions. 
# Usage = "QTG_Finder_predict.py -gl QTL_gene_list -sp species_abbreviation"
# QTL_gene_list: this is the list of QTL genes to be ranked. QTLs should be seperated by '//'. Structure for each QTL: start with '//' and QTL name, then append a gene list. See 'SSQ_batch_QTL_genes.csv' for a example  
# species_abbreviation: "AT" for Arabidopsis; "OS" for rice
# Usage example: python QTG_Finder_predict.py -gl SSQ_batch_QTL_genes.csv -sp 'AT'

import numpy as np
import random
import pandas as pd
from sklearn import model_selection  
from sklearn import metrics
from sklearn import ensemble
import time
import itertools as it
import argparse
import pickle

parser = argparse.ArgumentParser(description='QTL_Finder_v1.0')
requiredArg = parser.add_argument_group('required arguments')
requiredArg.add_argument('-gl',dest='gl',required=True, help='Input QTL gene list: this is the list of QTL genes to be ranked. See SSQ_batch_QTL_genes.csv for a example')
requiredArg.add_argument('-sp',dest='sp',required=True,choices=['AT', 'OS'], help='Species: use "AT" for Arabidopsis; "OS" for rice')
args = parser.parse_args()

def train_qtg(df, Validation_set_i):

    prediction_list = (len(Validation_set_i))*[0]
    validation_feature_all = Validation_set_i.drop(['class'], axis=1)  
    neg_inter= pickle.load(pik_f) #  number of iter models
    print(neg_inter)      
    for i in range(0,neg_inter ): 
        single_model=pickle.load(pik_f) # unpack individual models
        validation_pred = single_model.predict_proba(validation_feature_all)[:,1] # predictions based on models
        prediction_list+= validation_pred # combine prediction results
    prediction_list_all=np.array(prediction_list) 
    Ranks=prediction_list_all.argsort()[::-1].argsort() # rank genes based on times of being predicted as a causal gene
    prediction_list_freq= [l/(neg_inter) for l in prediction_list] # calculate freqency of being predicted as a causal gene
    rank_list_df= pd.DataFrame( {'ID':Validation_set_ID_uni, 'freq': prediction_list_freq }) 
    for i in gene_ex: # label excluded IDs with a freqency of 0
        rank_list_df= rank_list_df.append({'ID': i, 'freq':0 }, ignore_index=True) 

    rank_list_df['Rank'] = rank_list_df['freq'].rank(ascending=0,method='average') # For genes with the same frequency, the average rank is used for all of them
    rank_list_df_sorted=rank_list_df.sort_values(by=['Rank'],ascending =True) 
    rank_list_df_sorted=rank_list_df_sorted.reset_index(drop= True)

    with open('QTL_gene_rank.csv','a') as indenti_f: # output gene rank and frequency 
        indenti_f.write ('//'+'\n'+QTL_name+'\n')
        indenti_f.write ('ID'+','+'Rank_in_a_QTL'+','+'Frequency'+'\n')
        for i in range(len(rank_list_df_sorted)):            
            indenti_f.write (rank_list_df_sorted['ID'][i]+','+str(rank_list_df_sorted['Rank'][i])+','+ str(rank_list_df_sorted['freq'][i]) +'\n')

#### check known QTL gene on QTLs, remove, reducant ID, run
if __name__ == '__main__':
    start_time = time.time()
    with open(args.gl,'r') as f:  
        for key,group in it.groupby(f,lambda line: line.startswith('//') ): # QTLs are seperated by // in input file. 
            if not key:
                with open(args.sp+"_model.dat", "rb") as pik_f:
                    df = pickle.load(pik_f) # feature list with known causal genes used for training
                    group = list(group) # each group is a QTL
                    group=[i.strip('\n') for i in group]
                    QTL_name=group[0]
                    print('QTL: '+QTL_name)
                    genes_in_QTL=group[1:] # gene list
                    original_length= len(genes_in_QTL) # total number of input IDs in a QTL 
                    Validation_set=pd.DataFrame()   
                    for i in range(len(genes_in_QTL)):   
                        Validation_set=Validation_set.append (df[df.ID==genes_in_QTL[i]]) # take out IDs that can not be found in the feature list 

                    df = df.drop(['ID'], axis=1) 
                    ind_for_exclusion=[]
                    for t in range(len(Validation_set.index)): # exclude known causal genes that have been used for training 
                        if Validation_set['class'][Validation_set.index[t]]==1:
                            ind_for_exclusion.append(t)
                            print ('Known QTL causal gene excluded: ')
                            print(Validation_set['ID'][Validation_set.index[t]])  
                    Validation_set=Validation_set.drop(Validation_set.index[ind_for_exclusion])  # remove known causal genes from input gene list
                    Validation_set_ID_uni=list(Validation_set.ID) # unique ID, remove redandancy  
                    print('Number of genes: '+str(original_length)  )  # total number of input IDs in a QTL 
                    gene_ex=(set(genes_in_QTL)-set(Validation_set_ID_uni)) # exclude IDs that could not be found in feature list
                    Validation_set=  Validation_set.drop(['ID'], axis=1)  
                    Validation_set_i=Validation_set.reset_index(drop=True) # list of genes to be ranked
                    #train_set = df[df['class']==1]# train_set
                    train_qtg(df, Validation_set_i) # execute 
    print("--- %s seconds ---" % round((time.time() - start_time),2) ) 

