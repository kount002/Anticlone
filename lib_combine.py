# -*- coding: utf-8 -*-
"""
Created on Tue Jul 19 13:34:33 2016
Allows to work with clone counts from different pikles
Need to supply path to two pickled pandas libraries
Usage: lib_combine.py ~/anal/path2 ~/anal/path2
"""

import sys
import pickle
import pandas as pd
import numpy as np
import matplotlib.pylab as plt

def open_lib(path):
    try: 
        with open(path1, 'rb') as f:
            lib=pickle.load(f)
    except:
        print('Cannot open', path1, '\n check path')
    return(lib)



##############main##############
path1=sys.argv[1]
path2=sys.argv[2]

#parameters 
path1='~/anal/human/hiseq_trio/clone_count/clone_count_df.pkl'
path2='~/anal/human/phm_diversity/custom_prime/clone_count/clone_count_df.pkl'
keep=['A16000']


#join tables
df1=open_lib(path1)
df2=open_lib(path2)
print('Columns from right table:', df1.columns)
print('Columns from left table:', df2.columns)
df=df1.join(df2[keep], how='outer', lsuffix='l')




