# -*- coding: utf-8 -*-
"""
Created on Tue Jul 19 13:34:33 2016
Allows to work with clone counts from different pikles
Need to supply path to two pickled pandas libraries
Usage: lib_combine.py ~/anal/path2 ~/anal/path2
"""

import matplotlib
matplotlib.use('Agg')
import sys
import pickle
import pandas as pd
import numpy as np
import matplotlib.pylab as plt
from sklearn import linear_model
import re

def open_lib(path): #loads pickled df into memory
    try: 
        with open(path, 'rb') as f:
            lib=pickle.load(f)
        print('Read', path)
        return(lib)
    except:
        print('Cannot open', path, '\n check path')



def norm_varr(df):
# clone_count.py outputs raw clone counts. Function normalizes count as: count*ave(sum of sums)/specific sums
#may need to insert NaN back in place of 0. See what are your needs will be
#log transfom in the last expression    
    #very fast compared to pandas version

    lc=list(df.columns)
    lcs=[x for x in lc if x!='Annotation']
    df.fillna(0, inplace=True)
    aar=df[lcs].values
    ssums=aar.sum(axis=0) #sum based on column
    amp=ssums.mean()
    aar=(aar*amp)/ssums
    df[lcs]=pd.DataFrame(aar, index=df.index)
    df.replace(0, np.nan, inplace=True) #replaces all zeros for NaN
    df[lcs]=np.log10(df[lcs]) #log-transfrom data
    
    return(df)

def anal_prep(df):
    #prepares the dataframe for expression analysis, log transform it
    #changes all values <10 to 10
    #removes clones with all log(10)
    lc=list(df.columns)
    lcs=[x for x in lc if x!='Annotation']
    dfx=df.fillna(0)
    arr=np.array(dfx[lcs])
    arr[arr<1]=1
    dfx[lcs]=pd.DataFrame(arr, index=dfx.index)
    dfx.replace(to_replace=1, value=np.nan, inplace=True) # prepares for removal of clones with all basal expression (1)
    dfx.dropna(subset=lcs, how='all', inplace=True)
    dfx.replace(to_replace=np.nan, value=1, inplace=True)    

    return(dfx)


##############main##############
#path1=sys.argv[1]
#path2=sys.argv[2]

#parameters 
path2='/home/kount002/anal/human/miphage/clone_count/clone_count_df.pkl'
path1='/home/kount002/anal/human/phm_diversity/custom_prime/clone_count/clone_count_df.pkl'
keep=['A16000', 'Annotation']


#join tables from two df pickles
df1=open_lib(path1)
df2=open_lib(path2)
print('Columns from right table:', df1.columns)
print('Columns from left table:', df2.columns)

df=df1.join(df2[keep], how='outer', lsuffix='l')
df.dropna(how='all', axis=1, inplace=True) #drop empty rows
df['Annotation'].fillna(df.Annotationl, inplace=True) #collaple Annotations into one column
del(df['Annotationl'])
df.replace(to_replace=np.nan, value=0, inplace=True) #prep data for regression and lambda

mod=(lambda x: True if re.search(r'EK[0-9]{3}(?!N)', x) else False)
xl=[x for x in df.columns if mod(x)]


#normalizes values and log transforms them before regression
df=norm_varr(df)
df=anal_prep(df)


#prep data for multiple regression
dfy=df.A16000
dfx=df[xl]

model=linear_model.LinearRegression(fit_intercept=False)
model.fit(dfx, dfy)

coef=model.coef_
print('Coefficients for ', dfx.columns, 'are', coef)

#plot scatter for predicted model

py=model.predict(df.iloc[:,:3])

with open('560_561_570LR.pkl','wb') as f:
    pickle.dump(py, f)
 
plt.clf()
plt.scatter(dfy, py)
plt.savefig('test.png')

