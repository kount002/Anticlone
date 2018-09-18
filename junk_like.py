# -*- coding: utf-8 -*-
"""
Created on Fri Sep  8 13:28:29 2017

@author: Evgueni
"""
'''
test identification of the genes that while have lower count
but behaive like background genes
'''

################## Param ############
junklist=['COL4A2', 'NCOA6', 'MYCBP2', 'MAP1B', 'MAP1A', \
'TCEB3', 'PTN', 'USP10', 'ANKRD11', 'MTRNR2L2', 'ACTB', 'ACTG1'
'POLR3C','BCLAF1','PHF3', 'NCAN']

rutlist=['KRut120', 'KRut220', 'KAbMix20', 'ERut10', 'ERut11', 'ERut12', \
'ANA_1', 'ANA_17']

noslist=['KNOS120', 'KNOS20', 'KNOS220', 'E153']

file ='../ANA/Combine_nextseq/norm_by_PNTextraplus.csv'
file='../Tri_stage/hts_count_merged.txt'
file='../ANA/Combine_nextseq/next_hts_count_merged_recount.txt'
file='../170930_MSseq/Next5_norm_contgene_central.csv'
file='../170930_MSseq/Next5_quantile.csv'

################## Imports ##########
import pandas as pd
import numpy as np
import seaborn as sns
import explore_graph as exg

################## Func #############
def open_file(file): 
    ''' open file checks and removes Annotation column,
    not_aligned reads, all zeros'''

    if file.endswith('txt'):
        df=pd.read_csv(file, header=0, index_col='gene', sep=' ')
    elif file.endswith('csv'):
        df=pd.read_csv(file, index_col='gene')
    try:
        del(df['Annotation'])
    except:
        pass
    if str(df.iloc[:-1,:].index).endswith('not_unique'):
        df=df.iloc[:-5,:] #remove not aligned reads
        df=df[df.iloc[:,:].sum(axis=1)>0] #cuts genes with zero reads in all samples

    return(df)    

def low_cut(df, noslist, cut=1): 
    '''removes genes that have counts bellow provided parameter in 
    ALL samples 
    Needs a list of NOS samples to exclude from check'''
    
    df=df.loc[(df[df.columns.difference(noslist)]>cut).any(axis=1),:] #remove genes with all 0s or adjust low counts
    return(df)

def rut_cut(df, rutlist, noslist, deviation=2):
    ''' removes genes that have counts below provided number of deviations
    from mean of rutaxin controls
    For mean calculation zeros are excluded
    Needs a list of rutaxin samples
    Needs a list of NOS samples to exclude'''
    df['rut_mean']=df[rutlist][df[rutlist]>0].mean(axis=1) #excludes zero values
    df['rut_dev']=df[rutlist][df[rutlist]>0].std(axis=1)    #excludes zero values
    df['rut_cutoff']=df['rut_mean']+df['rut_dev']*deviation
    df.fillna(0, inplace=True)
    uselist=[x for x in list(df.columns) if x not in noslist] #restrict columns to human samples and removes control samples
    arr=df[uselist].values #converts to np.array
    boolseries=arr.T>arr[:,-1] #compares each column with the last column; produces a boolean array
    boolseries=np.any(arr.T>arr[:,-1], axis=0) # flatten bool array to one dimention must have one True to be True
    df=df.iloc[boolseries,:] 
    df=df.drop(['rut_cutoff','rut_mean','rut_dev'], axis=1)
    return(df)
    

    
################## Main #############


#dfcorr=df.T.corr()
#dfsel=dfcorr.ix['ACTB']
#sns.clustermap(dfcorr.ix[junklist, junklist].fillna(0))
#dfn=exg.pre_norm(df)
 

df=open_file(file)
df=low_cut(df, noslist)
dfn=rut_cut(df, rutlist, noslist, deviation=3)


