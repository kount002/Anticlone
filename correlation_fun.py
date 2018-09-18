# -*- coding: utf-8 -*-
"""
Created on Tue Oct  3 16:03:09 2017

@author: Evgueni
"""
################ Param ##############

contgene=['COL4A2', 'NCOA6', 'MYCBP2', 'MAP1B', 'MAP1A', \
'TCEB3', 'PTN', 'USP10', 'ANKRD11', 'MTRNR2L2', 'ACTB', 'ACTG1',\
'POLR3C', 'PHF3', 'NCAN', 'NCOA6', 'MTRNR2L8', \
'PPIA']

extgene=['MALAT1',  'ANKRD1']
mincorr=-0.3
maxcorr= 0.3

file='../170930_MSseq/Junknorm_comb.csv'

############## Imports #############
import pandas as pd
import explore_graph as exg

############## Funct  ##############


############# Main ##############
df=pd.read_csv(file, header=0)
#df=df.iloc[:-5,:] #drop no feature junk at the bottom
df.set_index('gene', drop=True, inplace=True)
#lc=df.columns
#lc=[x.strip('sample_') for x in lc]
#df.columns=lc
dflight=exg.corr_cull(df, contgene, extgene, mincorr, maxcorr)

