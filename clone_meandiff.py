# -*- coding: utf-8 -*-
"""
Created on Mon Jul 25 11:04:37 2016
USAGE: clone_meandiff.py master.pkl
where master.pkl is not normalized pickle

Determines means of sample groups than subtract from them 
max(mean(control1),mean(control2))
Sort items by max expression or max difference between
2 x mean(interest)-(mean(gr1)+mean(gr2))

Uses graphs that are sample specific, 
"""
################### Params ###############
inputs='master.pkl' #path
output='diff_out_master.csv' #path
#fname='p_value_hist.png' # path to figure

method='upper85' #normalalizaiton method (upper80, tc, RLE80, med, max)
tresh=1 #normalization based on minimum absolute count (provide minimum count)
meanfilter=1.3 #normalization based on mean expression (folds of the tresh) (removes genes wt mean expression less then value)


foldover=1 #fold over max(control)

groups={
    'NMO_interest':['K20120','K20320','K20420','K20520','K20620'],
    'Healthy_compare':['K10120','K10220','K10320','K10420','K10520'],
    'NOS_control':['KNOS120','KNOS220','KNOS20'],
    'RUT_control':['KRut120','KRut220','KAbMix20']
        }

#'K10820','K10920','K11020','K11120','K11220'],
#['K20120','K20320','K20420','K20520','K20620'
#'K10120','K10220','K10320','K10420','K10520'
#'KRut120','KRut220','KAbMix20'


############ imports #####################
import pickle
import pandas as pd
import numpy as np
import explore_graph as exg



############# functions ###################
def content(c):
    ''' prints contents of the c (df) '''
    for i in c:
        print(i)

############# main ########################
df=exg.open_lib(inputs)
print ('Opened sample set {0}, containing:'.format(inputs) )
cont=df.columns
content(cont)

#check if named samples are present in the library
samples=[y for x in groups.values() for y in x] #flatten groups dict 
for i in samples:
    if i not in cont:
        print(i, 'Not present in the sample library, check PARAM section')

#asks if the file works with gene list or fragment list and reduces bins to single position for fragment file        
if 'gene' not in list(cont):
    print('To collapse change line 68: Collapsing fragment bins to single-end bins')
    #df=exg.single_end(df)

    
##clean up annotation
df=exg.annot_clean(df)

#remove items with no feature in alingment position
print('Array shape for all clones', df.shape)
df=df[~df['Annotation'].str.startswith("__")]
print('Array shape after unnotated clones removed', df.shape)

#keep only that are named in parameter section
keeps=samples+['Annotation']
df=df[keeps]
#run normalizatiion routine
df=exg.norm_varr(df, method, tresh, meanfilter)

#create df with means of the groups
dfm=pd.DataFrame()
for k, i in groups.items():
    nname='mean_'+k
    dfm[nname]=df[i].mean(axis=1)
    dfm[i]=df[i]

dfm['Annotation']=df['Annotation']
cols=dfm.columns #put Annotation in front of all columns
cols=cols[-1:]+cols[:-1]
dfm=dfm[cols]

contr=['mean_'+x for x in groups.keys() if x.find('control')>=0]
dfm['max_control']=dfm[contr].max(axis=1)

#filter out all that less than max(controls)
interkeys=[x for x in dfm.columns if x.find('interest')>=0]
dfmb=dfm.loc[dfm[interkeys[0]]>foldover*dfm['max_control']]

dfmb=dfm # does not do control subtraction, only compares to Healthy
    
#filter out all that less than max(controls)
##clean up annotation
#dfmb=exg.annot_clean(dfmb)

#plots to check for normalization using unculled list
if len(dfm)>10000:
    dfms=dfm.sample(10000)
else:
    dfms=dfm

exg.MA_plot(dfms[interkeys[0]],dfms['mean_NOS_control'], 'mafile')
exg.expression_plot(dfms[interkeys[0]],dfms['max_control'], 'exfile')

if len(df)>10000:
    dfs=df.sample(10000)
else:
    dfs=df
exg.norm_plot(dfs)

#select what is over the Healthy control
comparekeys=[x for x in dfm.columns if x.find('compare')>=0]

#build two compare code
for i in comparekeys:
    nname='diff_'+i
    #dfmb.is_copy=False
    dfmb[nname]=dfmb[interkeys[0]]-dfmb[i]
#dfmb['diff_N']=dfmb[interkeys[0]]-dfmb[comparekeys[0]]
comparename=[x for x in dfmb.columns if x.find('diff')>=0]
dfmb.sort_values(comparename[0], axis=0, ascending=False, inplace=True)
dfmb.to_csv(output)

#testing
#a=dfm.loc[dfm['mean_NMO_interest']>dfm['mean_RUT_control']]
#a=dfm.loc[dfm['mean_NMO_interest']>dfm['max_control']]
#print(a.shape)
    
#print(dfm.shape)
print('Done')
