# -*- coding: utf-8 -*-
"""
Created on Fri Jul 22 16:16:00 2016

@author: Evgueni
"""

'''
hts_count_merge.txt
analysis of gene counts
'''
#########Parm#########
file='hts_count_merged.txt'
file='hts_comb_wt_mi_div_tr.txt'
file='hts_hybrid.txt'

######################

import pandas as pd
import pickle
import explore_graph as exg
 
df=pd.read_csv(file, sep=' ', header=0)
df=df.iloc[:-5,:] #drop no feature junk
df.set_index('gene', drop=True, inplace=True)
lc=df.columns
lc=[x.strip('sample_') for x in lc]
df.columns=lc
df['Annotation']=df.index #mods to accomodate gene lists

df=exg.norm_varr(df, 'RLE90', tresh=1, meanfilter=1.0) #comment if run anova to avoid double normalization
exg.norm_plot(df)
exg.expression_plot(df.iloc[:,19], df.iloc[:,12], 'second')
#exg.MA_plot(df.iloc[:,18], df.iloc[:,16], 'first')

df.to_csv('hts_hybrid.csv')
with open('mastera.pkl', 'wb') as f:
    pickle.dump(df, f)
    
#a=dfx.loc[(dfx>1.5).any(axis=1),:] #extract rows with values all>1.5
    

