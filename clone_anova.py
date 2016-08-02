# -*- codiing: utf-8 -*-
"""
Created on Mon Jul 25 11:04:37 2016
USAGE: clone_stat.py master.pkl
Runs ANOVA and Turkey tests on ea gene/clone
and returns clones that are significantly different 

"""
################### Params ###############
inputs='master.pkl' #path
save='stat_out_master.csv' #path
#check list of analyses?
kruskal=0 #if set to 1 will use non-parametric anova
method='med' #normalization: 'med' for median, 'tc' for total count
groups={
    'gr1':['K20120','K20320','K20420','K20520','K20620'],
    'gr2':['K10820','K10920','K11020','K11120','K11220'],
    'gr3':['K10120','K10220','K10320','K10420','K10520']
        }

#'K10820','K10920','K11020','K11120','K11220'],
#['K20120','K20320','K20420','K20520','K20620'
#'K10120','K10220','K10320','K10420','K10520'


############ imports #####################
import sys 
import pandas as pd
import numpy as np
import scipy.stats as stats
import explore_graph as exg
import pickle
#import statmodels.formula.api as ols


############### functions ################

def open_lib(path): #loads pickled df into memory
    try: 
        with open(path, 'rb') as f:
            lib=pickle.load(f)
        print('Read', path)
        return(lib)
    except:
        print('Cannot open', path, '\n check path')
        sys.exit(2)
        
def convert_gr(df, groups):
    '''converts the provided dictionary with columns to a dic
    with numerical indices for the corresponding columns
    needed for use with itertuples
    '''
    groupn={}
    for cols in groups:
        tm=[]
        for col in groups.get(cols):
            idx=df.columns.get_loc(col)
            tm.append(idx)
        groupn[cols]=tm #dic with corresponsind indeces   
    return(groupn)    

def get_values(dftup, groupn):
    ''' assign expression values to a list of lists for use as *args '''
    args=[]
    for i in groupn.keys():
        ls=[dftup[x] for x in groupn[i]]
        args.append(ls)
    return(args)
        
#################### main ###################
def main():        
    print('Using parametric/non-parametric Anova "0/1:"', kruskal)
            
    #load library
    df=open_lib(inputs)  #load library
    df=exg.norm_varr(df, 'med') #normalize and log transform
    df.fillna(0, inplace=True)
    
    #convert groups ID to indeces
    groupn=convert_gr(df, groups)
    
    #Anova one way test 
    Flist=[]
    plist=[]
    for i in df.itertuples():
        args=get_values(i[1:], groupn)
        if not kruskal:
            F,p=stats.f_oneway(*args) #parametric Anova
        else:
            try: F,p=stats.mstats.kruskalwallis(*args) #non-parametric
            except: pass
        Flist.append(F)
        plist.append(p)    
        
    arrt=np.array([Flist, plist]).T
    ds=pd.DataFrame(arrt, index=df.index, columns=['F_statistic', 'p_value'])
    df=df.join(ds)
    #extract only significat genes
    dfsig=df[df['p_value']<0.05].sort_values('p_value')[df['p_value']<0.05]
    
    with open(save, 'wb') as f:
        dfsig.to_csv(save)


if __name__ =='__main__':
    main()

        
