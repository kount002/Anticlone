# -*- coding: utf-8 -*-
"""
Created on Tue Jul 26 15:41:48 2016
USAGE: lib_combine [*]
if any agrument present the scipt will just output listing of columns
though it will create a dummy master.pkl

input data (pickled dataframes and samples) are organized as dictionary:
path:[samples]
empty sample == all

"""
import pandas as pd
import numpy as np
import pickle
import sys


################ Parameters #############
#empty will use all columns
inputs={
    
    'clone_count/clone_count_df.pkl':
    [],
    '../phm_diversity/custom_prime/560_561_570LR.pkl':
    [],
    '../miphage/clone_count/clone_count_df.pkl':
    ['G20RUT', 'E20201'] ,
    '../tri_stage/clone_count/clone_count_df.pkl':
    ['20120','AbMix20', 'NOS20']
 
    
    }
save='master.pkl' #name of file to save to

################# functions ############

def open_lib(path): #loads pickled df into memory
    try: 
        with open(path, 'rb') as f:
            lib=pickle.load(f)
        print('Read', path)
        return(lib)
    except:
        print('Cannot open', path, '\n check path')
        sys.exit(2)


################# main ##################


dfmast=pd.DataFrame()  #make master df from 
dfmast['Annotation']=np.nan
for i in inputs.keys():
    dft=open_lib(i)
    print('Columns from right table:', dft.columns)
    if len(sys.argv)>1: #test if user asked for column names only
        continue
    keep=inputs.get(i)
    if not keep:
        keep=dft.columns
    else:
        keep.append('Annotation')
    dfmast=dfmast.join(dft[keep], how='outer', lsuffix='l')
    dfmast['Annotation'].fillna(dfmast['Annotationl'], inplace=True) #collaple Annotations into one column
    del(dfmast['Annotationl'])
    
print('\n Master df was build. Shape:{0}; Columns {1}'.format(dfmast.shape, dfmast.columns))
with open(save, 'wb') as f:
    pickle.dump(dfmast, f)
print('\n Saved as', save)




