# -*- coding: utf-8 -*-
"""
Created on Tue Jul 19 13:34:33 2016
Allows to work with clone counts from different pikles
Need to supply path to two pickled pandas libraries
Usage: lib_combine.py ~/anal/path2 ~/anal/path2
Requires explore_graph.py in the same folder to work
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


def norm_f(df):
    lc=list(df.columns)
    lcs=[x for x in lc if x!='Annotation']
    df.replace(0, np.nan, inplace=True) #replaces all zeros for NaN
    df[lcs]=np.log10(df[lcs]) #log-transfrom data
    return(df)


##############main##############
#path1=sys.argv[1]
#path2=sys.argv[2]

#parameters 
path2='/home/kount002/anal/human/miphage/clone_count/clone_count_df.pkl'
path1='/home/kount002/anal/human/phm_diversity/custom_prime/clone_count/clone_count_df.pkl'
keep1=[] #empty will use all columns
keep2=['A16000', 'Annotation'] #empty will use all columns


#join tables from two df pickles
df1=open_lib(path1)
df2=open_lib(path2)
print('Columns from right table:', df1.columns)
print('Columns from left table:', df2.columns)

if not keep1:
    keep1=df1.columns
if not keep2:
    keep2=df2.columns


df=df1[keep1].join(df2[keep2], how='outer', lsuffix='l')
df['Annotation'].fillna(df.Annotationl, inplace=True) #collaple Annotations into one column
del(df['Annotationl'])
mod=(lambda x: True if re.search(r'EK[0-9]{3}(?!N)', x) else False)
xl=[x for x in df.columns if mod(x)]

dfr=df.copy(deep=True)
#remove low values and log transforms them before regression use filtered values for model and real for prediction
df=np.log10(df[xl+['A16000']])
df.replace(np.nan, 1, inplace=True)
dfm=df.loc[(df.iloc[:,:4]>1).any(axis=1),:] #extract rows with values over 1


#prep data for multiple regression
dfy=dfm.A16000
dfx=dfm[xl]

model=linear_model.LinearRegression(fit_intercept=False)
model.fit(dfx, dfy)

coef=model.coef_
print('Coefficients for ', dfx.columns, 'are', coef)

dfpy=dfr.iloc[:,:3].pow(coef, axis=1).product(axis=1) #formula to calculate predicted value based on coef
#sum((expression count)**coef)
dfpy=pd.DataFrame(dfpy)
dfpy.columns=['Regress']
dfpy['Summed']=dfr.iloc[:,:3].sum(axis=1)
dfpy=dfpy.join(dfr[['A16000', 'Annotation']], how='inner')

with open('560_561_570LR.pkl','wb') as f:
    pickle.dump(pd.DataFrame(dfpy), f)
#plot scatter for predicted model

dfpy.replace(np.nan, 0, inplace=True)
#py=model.predict(df.iloc[:,:3]) #adjust fo reduced counts after small coef


plt.figure(figsize=(6,6))
plt.scatter(np.log10(dfr['A16000']), np.log10(dfpy['Summed']))
plt.savefig('Scatter_test.png')
 
