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

contgene=['COL4A2', 'NCOA6', 'MYCBP2', 'MAP1B', 'MAP1A', \
'TCEB3', 'PTN', 'USP10', 'ANKRD11', 'MTRNR2L2', 'ACTB', 'ACTG1',\
'POLR3C', 'PHF3', 'NCAN', 'NCOA6', 'MTRNR2L8', \
'PPIA']

file='hts_count_merged.txt'
file='hts_count_merged_recount.txt'
file='../ANA/Combine_nextseq/comb_hts.txt'
file='../170930_MSseq/Raw_comb.csv'


######################

def pre_norm (file): #place holder, the function is impemented in explore graph

    '''Runs prenormalization for a set of genes and helps to select
    genes to use as normalization set. Those gene a junk, background
    genes that are assumed to be equivalent between the samples.
    Build on an assumptions that junk is abandunt, high count, and
    correlates in intensity between samples. 
    Function outputs a df with gene list and clustermap figure.
    Function takes a raw counts file, that has no Annotation column.
    '''
    
    import numpy as np
    import seaborn as sns
    
    dfraw=pd.read_csv(file)
    dfraw.replace(0, np.nan, inplace=True)
    dflog=np.log10(dfraw.iloc[:,:])
    dflog.index=dfraw['gene']
    dflog.fillna(0, inplace=True)
    dflogs=dflog[dflog.iloc[:,:].min(axis=1)>1.6]
    sns.clustermap(dflogs.T.corr())
    df.to_csv('norm_by_junt.csv')
    return(df)
    

def corr_cull (df, contgene, extgene, mincorr, maxcorr):
    ''' removes genes that correlate with PTN and MALAT1 junk.
    takes a dataframe, contgene - list of genes used for normlization; 
    extgene -  list of genes that are non-selected junk like MALAT1;
    mincorr, maxcorr - min and max correlation cutoffs
    '''
    glist=contgene+extgene #combine list of genes to correlate with
    try:
        del(df['Annotation'])
    except:
        pass
    dfcorr=df.T.corr()
    dflight=df
    for i in glist:
        dflight=dflight[ (dfcorr[i]>=mincorr) & (dfcorr[i]<=maxcorr)]
    return(dflight)
    
    
######################

import pandas as pd
import pickle
import explore_graph as exg
 
df=pd.read_csv(file, header=0)
df=df.iloc[:-5,:] #drop no feature junk at the bottom
df.set_index('gene', drop=True, inplace=True)
lc=df.columns
lc=[x.strip('sample_') for x in lc]
df.columns=lc
df['Annotation']=df.index #mods to accomodate gene lists



df=exg.norm_varr(df, 'contgeneRLE75', contgene, tresh=10, meanfilter=1.1) #comment if run anova to avoid double normalization
exg.norm_plot(df)
#exg.expression_plot(df.iloc[:,19], df.iloc[:,12], 'second')
#exg.MA_plot(df.iloc[:,18], df.iloc[:,16], 'first')

#with open('master.pkl', 'wb') as f:
#    pickle.dump(df, f)
    
#a=dfx.loc[(dfx>1.5).any(axis=1),:] #extract rows with values all>1.5
    
df.to_csv('../170930_MSseq/Junknorm_comb.csv')
