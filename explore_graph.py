# -*- coding: utf-8 -*-
"""
Created on Thu Jun 30 20:15:09 2016

@author: Evgueni
"""

import sys
import pandas as pd
import numpy as np
import pickle
import matplotlib
matplotlib.use('Agg')  #need this to turn off figure display
import matplotlib.pyplot as plt




def norm_df(df):
# clone_count.py outputs raw clone counts. Function normalizes count as: count*ave(sum of sums)/specific sums
    lcs=[x for x in lc if x!='Annotation']
    for col in lcs:
        print('Processing...', col )
        df[col]=df[col].apply(lambda x: x/df[col].sum()*df[lcs].sum().mean())
    return(df)
    
def norm_varr(df):
    #very fast compared to pandas version
    lcs=[x for x in lc if x!='Annotation']
    df.fillna(0, inplace=True)
    aar=df[lcs].values
    ssums=aar.sum(axis=0) #sum based on column
    amp=ssums.mean()
    aar=(aar*amp)/ssums
    df[lcs]=pd.DataFrame(aar, index=df.index)
    
    return(df)


########################## MAIN  ###########################

if len(sys.argv)<4: #check for agrument presence
    print('Usage: graphing.py pickle_df.pkl base_sample \
    compare_sample [compare_sample2 ...]')
    sys.exit(2)

try:    #check that pickle exists
    with open(sys.argv[1], 'rb') as handle:
        df=pickle.load(handle)
except:
    print('Cannot open pickle file. Run clone_count or verify the path.')
    sys.exit(2)


lc=list(df.columns)

print('The DataFrame contains following columns:', lc)

#dfdf=norm_df(df)
dfar=norm_varr(df)
