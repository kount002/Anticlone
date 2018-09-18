# -*- coding: utf-8 -*-
"""
Created on Mon Jul 25 11:04:37 2016
USAGE: clone_stat.py master.pkl
Runs ANOVA and Turkey tests on ea gene/clone
and returns clones that are significantly different 

"""
################### Params ###############
inputs='hts.pkl' #path
inputs='../170930_MSseq/Junknorm_comb.csv'
save='stat_out_master_HN.csv' #path
fname='p_value_hist.png' # path to figure
method='upper80' #normalalizaiton method (upper, tc, med, max)
#check list of analyses?
kruskal=1 #if set to 1 will use non-parametric anova
welch=0 #(use with kruskal=0); if set to 1 will use Welch if '0' will ignore
groups={
    'gr1':['K10520','K10920','K11020','K11120','K11220', 'E107', 
           'E113', 'E114', 'E115', 'E116', 'E117'],
    'gr2':['K10120','K10220','K10320','K10420','K10820', 'E150',
           'E151', 'E152', 'E153', 'E000','EInLib']
        }

#'K10820','K10920','K11020','K11120','K11220'],
#'K10120','K10220','K10320','K10420','K10520'
#'KRut120','KRut220','KAbMix20'
#'K20120','K20320','K20420','K20520','K20620', 'E202', 
#    'E207', 'E208', 'E209', 'E210', 'E211'
# 'K10120','K10220','K10320','K10420','K10820', 'E107', 
#    'E150', 'E151', 'E152', 'E153', 'E000','EInLib'


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
    print('Using parametric/non-parametric Anova "0/1"', kruskal)
            
    #load library
    if inputs.endswith('pkl'):
        print( 'Using pickle input not, will normalize')
        df=open_lib(inputs)  #load library
        df=exg.norm_varr(df, method) #normalize and log transform
        df.fillna(0, inplace=True)

    else:
        print('Using csv input')        
        df=pd.read_csv(inputs, header=0)    
    #convert groups ID to indeces
    groupn=convert_gr(df, groups)
    
    #Anova one way test 
    Flist=[]
    plist=[]
    for i in df.itertuples():
        args=get_values(i[1:], groupn)
        if not kruskal:
            if welch:
                F,p=stats.ttest_ind(args[0], args[1], equal_var=False)
            else:
                F,p=stats.f_oneway(*args) #parametric Anova            
        else:
            F,p=np.nan, np.nan
            try: F,p=stats.mstats.kruskalwallis(*args) #non-parametric
            except: pass
        
        Flist.append(F)
        plist.append(p)    
        
    arrt=np.array([Flist, plist]).T
    ds=pd.DataFrame(arrt, index=df.index, columns=['F_statistic', 'p_value'])
    df=df.join(ds)
    #extract only significat genes
    dfsig=df[df['p_value']<0.05].sort_values('p_value')[df['p_value']<0.05]
    fig=df.iloc[:,-1].hist(bins=40)
    fig=fig.get_figure()
    fig.savefig(fname)
    
    dfsig.to_csv(save)


if __name__ =='__main__':
    main()
