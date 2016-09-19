# -*- coding: utf-8 -*-
"""
Created on Wed Jul 27 18:54:20 2016
USAGE: clone_substract.py master.pkl
substracts expression level of negative controls from test samples
takes log-transfromed data, so works with fold increase in expression


"""
#################### Param ######################
inputs='master.pkl' #path to the dataframe
save='fold_out_master.csv' #path result file
savedf='subtract_master.pkl' #path to dataframe that has been reduced

method='tc' #choose normalization method between 'tc' and 'med'
#negative grOUPS to subtract
groups={
    'gr1':['KNOS120','KNOS220','KNOS20', 'H20NOS', 'NOS20' ],
    'gr2':['KRut120','KRut220', 'G20RUT'],
    'gr3':['AbMix20', 'KAbMix20'],
    'gr4':['Regress', 'Summed', 'A16000'],
    'gr5':['K10120','K10220','K10320','K10420','K10520']
    
        }
#positive groups to focus:
#focusgr={
#    'fgr1':[],
#    'fgr2':[],
#    'fgr3':[]
#        }




#################### Imports ####################
import sys
import pickle 
import pandas as pd
import numpy as np
import clone_anova as csa
import explore_graph as exg

#################### Functions ##################
def subtract_vals(df, groups):
    ''' outputs a column with subtract values calculated as max of 
    means in groups
    negative values are replaced with 0'''
    dft=pd.DataFrame(index=df.index)
    for gr in groups:
        tmean=df[groups[gr]].mean(axis=1) #mean function is here (change to max for stringency?)
        dft[gr]=tmean
    dft=dft.max(axis=1)
    dft[dft<0]=0
    return(dft)

def sub_list(groups):
    ''' compiles columns names from the dictionary to a single list '''
    sublist=[]
    for gr in groups:
        sublist+=groups[gr]
    return(sublist)
    
    
################### Main ########################


df=csa.open_lib(inputs) #open lib
df=exg.norm_varr(df, method)    #normalize and log transform
#df.fillna(0, inplace=True)

dft=subtract_vals(df, groups) #determine subtract value

sublist=sub_list(groups) #convert dic in to list
sumlist=[x for x in df.columns if x not in sublist] #make list of samples to treat
sumlist=[x for x in sumlist if x!='Annotation']
dfs=pd.DataFrame(index=df.index)
for col in sumlist: #?? change to sumlist/df.columns for comprehensive
    name=col+'_fd'
    dfs[name]=df[col]-dft
foldlist=[x+'_fd' for x in sumlist]
dfs=dfs.loc[(dfs[foldlist]>0).any(axis=1),:]

#join with Annotation column
try:
    dfs.join(df['Annotation'], how='left')
except:
    pass
#stop here if want fold increase over background table, move to save

#Use dfs table as a mask to cull df
dfsidx=dfs.index
dfs=df.loc[dfsidx]

####
dfs.to_csv(save)
with open(savedf, 'wb') as f:
    pickle.dump(dfs, f)
print('Scipt completed, see "{0}" for results.'.format(save)) 

