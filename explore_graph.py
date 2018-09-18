''' script to test data using pandas dataframe

'''
import sys
import pandas as pd
import numpy as np
import pickle
import matplotlib
matplotlib.use('Agg')  #need this to turn off figure display
import seaborn as sns
import matplotlib.pyplot as plt

def open_lib(path): #loads pickled df into memory
    try: 
        with open(path, 'rb') as f:
            lib=pickle.load(f)
            print('Read', path)
        return(lib)
    except:
        print('Cannot open', path, '\n check path')
        sys.exit(2)

def annot_clean (df):
    '''cleans up annotation column by removing repeated items for each row'''
    def transf(annt): 
        annt=str(annt)
        for i in ['(', ')', '\'', '\"', ' ']:
            annt=annt.replace(i, '')
        annt=annt.split(',')
        annt=set(annt)
        annt=list(annt)
        annt=[x for x in annt if x!=''] 
        annt.sort()
        annt=''.join(annt[:1]) # remove if want to keep set of hits instead of single value
        #annt=', '.join(annt[0:]) #uncomment if wat to keep set of hits insted of single value 
        return(annt)
   
    df['Annotation']=df['Annotation'].apply(transf)
    return(df)
    
def single_end(df):
    '''colapses bins based on start and end position to
       only start position bins
       works with fragment library only'''

    loc=df.index
    df['Chr'], df['Pos1'], df['Pos2']=zip(*loc)
    del(df['Pos2'])
    df['Chr']=df['Chr']+', '+(df['Pos1'].astype('str'))
    del(df['Pos1'])
    lc=df.columns
    lcnum=[x for x in lc if x!='Annotation']
    dfnum=df[lcnum]
    dfobj=df[['Annotation', 'Chr']]
    dfnum=dfnum.groupby(['Chr']).sum()
    dfnum.index.name=None
    dfobj=dfobj.groupby(['Chr']).apply(lambda x: ', '.join({str(a) for a in x['Annotation']}))
    dfobj.index.name=None
    dfobj.columns=['Annotation']
    dfobj.name='Annotation'
    df=pd.concat([dfnum, dfobj], axis=1)
    return(df)

def pre_norm (dfraw):
    '''Runs prenormalization for a set of genes and helps to select
    genes to use as normalization set. Those gene a junk, background
    genes that are assumed to be equivalent between the samples.
    Build on an assumptions that junk is abandunt, high count, and
    correlates in intensity between samples. 
    Function outputs a df with gene list and clustermap figure.
    Function takes a raw counts df, that has no Annotation column
    and nor _ambious reads.
    '''
    
    #dfraw=pd.read_csv(file, header=0, sep=' ')
    dfraw.replace(0, np.nan, inplace=True)
    try:
        dfraw.index=dfraw['gene']        
        del(dfraw['gene'])  
    except:
        pass
    dflog=np.log10(dfraw.iloc[:,:])
    dflog.fillna(0, inplace=True)
    dflogs=dflog[dflog.iloc[:,:].min(axis=1)>1.1]
    fig=sns.clustermap(dflogs.T.corr())
    plt.setp(fig.ax_heatmap.get_yticklabels(), rotation=0)    
    dflogs.to_csv('norm_by_junt.csv')
    return(dflogs)

    
def norm_varr(df, method='tc', contgene=[], tresh=1, meanfilter=1.0):
    ''' clone_count.py outputs raw clone counts. 
    Takes values for normalization methods: 'tc' - total count or 'med' - median
    'tc': (default) Function normalizes count as: count*ave(sum of sums)/specific sums
    'med':  Function normalizes count as: count*ave(of medians)/specific median
    does ignores 0 counts for median calculation.
    'max': uses max value to normalize
    'upper': uses upper quartile normalization
    'contgeneRLE50': need contgene, and uses RLE50 for control gene set
    'quantile': uses a rank normalization
    Does not consider samples with 'NOS' in name when removes low count genes from analysis
    Can change filtration of low experssed genes, use lines at the end    
    
    
    #may need to insert NaN back in place of 0. See what are your needs will be
    #log transfom in the last expression    
    #very fast compared to pandas version '''

    def rle_norm(aar):
        '''subfunction that calculates rle for normalization.
        Created to minimize reuse of the code for overall and
        control gene list '''
        meansample=np.exp(np.mean(np.log(aar), axis=1)) #synthetic sample made of geo means
        meansample=meansample.reshape(-1,1)
        ratsample=aar/meansample #ratios between samples and synthetic
        supper=np.percentile(ratsample, perc, axis=0) #percientile value        
        amp=np.exp(np.mean(np.log(supper)))
        supper=supper.reshape(1,-1)
        return(amp, supper)

    lc=list(df.columns)
    lcs=[x for x in lc if x!='Annotation']
    lcsnos=[x for x in lcs if not x.upper().find('NOS')>=0] #exclude NOS samples (for reduction)
    #repls=[x for x in range(tresh)]  #filter out counts that <10. Do it here to add more weight to highly expressed genes
    #df.replace(repls, tresh, inplace=True) #filter out counts that <10
    #print(df.shape)    
    df.fillna(0, inplace=True)    
    df=df.loc[(df[lcs]>0).any(axis=1),:] #remove genes with all 0s or adjust low counts
    df.is_copy=False #turns off copy change warnin    
    aar=df[lcs].values
    #aar[aar<2]=0
    if method=='med': #DEPRECIATED, use upper50 instead
        print('Using median normalization.')
        aarm=np.ma.masked_array(aar, mask=aar<=0) #mask 0 counts
        smed=np.ma.median(aarm, axis=0) #median based on columns
        amp=smed.mean()
        aar=(aar*amp)/smed
    elif method=='max': #DEPRECIATED, use upper100 instead
        print('Using max value (max) normalization')
        smax=aar.max(axis=0) #max based on columns
        amp=smax.mean()
        aar=(aar*amp)/smax  
    elif method=='tc':
        print('Using total count (tc) normalization')
        ssums=aar.sum(axis=0) #sum based on column
        amp=ssums.mean()
        aar=(aar*amp)/ssums
    elif method.startswith('upper'):
        print('Using upper quartile normalization:', method)
        if method=='upper':
            perc=90
        else:
            perc=int(method.lstrip('upper'))
        aar=aar.astype(np.float) ##check performance Delete?
        aar[aar==0]=np.nan ##check performance and delete
        supper=np.nanpercentile(aar, perc, axis=0) #check, change nan
        print('Quartiles:', perc, supper)
        if 0 in supper:
            print('Encountered 0 value as percentile, use higher value for percentile')
        amp=np.exp(np.mean(np.log(supper))) #geometric mean centers on most frequent
        aar=np.nan_to_num(aar)
        aar=(aar*amp)/supper
    elif method.startswith('RLE'):
        print('Using RLE normalization:', method)
        if method=='RLE':
            perc=50
        else:
            perc=int(method.lstrip('RLE'))
        aar[aar==0]=1 #assume no nan left
        
        amp, supper=rle_norm(aar)
        '''
        meansample=np.exp(np.mean(np.log(aar), axis=1)) #synthetic sample made of geo means
        meansample=meansample.reshape(-1,1)
        ratsample=aar/meansample #ratios between samples and synthetic
        supper=np.percentile(ratsample, perc, axis=0) #percientile value        
        amp=np.exp(np.mean(np.log(supper)))
        supper=supper.reshape(1,-1)
        '''
        
        aar=(amp*aar)/supper
        
    elif method.startswith('contgeneRLE'):
        print('Using RLE of the control gene list (contgene var):', method)
        if method=='contgeneRLE':
            perc=50
        else:
            perc=int(method.lstrip('contgeneRLE'))
        #aar[aar==0]=1
        '''
        #intersect control gene list with master gene, in case of mismatches
        idxm=df[lcs].Index()
        idxc=pd.Index(contgene)
        idxm.intersection(idxc)
        '''
        
        contgene_d=pd.DataFrame(index=contgene) #empty df out of the control list
        
        #make a control gene count df
        contgene_df=pd.merge(contgene_d, df, how='inner', \
        left_index=True, right_on='Annotation')
        #contgene_df.index=contgene_df(['gene'])
        aarcont=contgene_df[lcs].values
        aarcont[aarcont==0]=1        
        
        amp,supper=rle_norm(aarcont)
        aar=(amp*aar)/supper        
        
    elif method.startswith('quantile'): 
        print('Using quantile rank normalization: replace nan', method)
        #df=10**df 
        df.replace(0, np.nan, inplace=True)
        rank_mean=df[lcs].stack().groupby(df[lcs].rank(method='first', ascending=False).stack().astype(int)).quantile(q=0.25)
        dfq=df[lcs].rank(method='first', ascending=False).stack().astype(int).map(rank_mean).unstack()
        dfq.fillna(0, inplace=True)
        aar=dfq.values
    else:    
        pass #for none normalization

    aar[aar<1]=1 #after norm some values go <1, and turn negative after log-transform
    df[lcs]=pd.DataFrame(aar, index=df.index)
    df=df.loc[(df[lcsnos]>=tresh).any(axis=1),:] #may adjust to vary filter stringency   
    #df[lcs]=df[lcs].where(df[lcs]>=tresh, other=tresh) #replaces values that are zero and <tresh with tresh value
    df.replace(0, np.nan, inplace=True) #replaces all zeros for NaN #no need if zeros replaced with tresh    
    df[lcs]=np.log10(df[lcs]) #log-transfrom data
    df.fillna(0, inplace=True) #reciprocates with (where) line 
    #df=df.loc[(df[lcs]>1).any(axis=1),:]
    #df=df.loc[(df[lcsnos]>1.0).any(axis=1),:] #may adjust value for a filter
    treshlog=np.log10(tresh*meanfilter) #placeholder for alternative filter on means
    #                                   calculates log of the minimum mean value that was set by user
    df=df.loc[(np.mean(df[lcsnos], axis=1)>treshlog),:] #alternative filter on mean across all samples    
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


def anal_prep(df):
    '''prepares the dataframe for expression analysis, log transform it
    changes all values <10 to 10
    removes clones with all log(10) '''
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

def norm_plot(df):
    ''' plot variance of all samples vs mean expression
        helps to determine which normalization works best for a set
        Creates two plots boxplot and scatter plot
    '''

    lc=list(df.columns)
    lcs=[x for x in lc if x!='Annotation']
    #lcsnos=[x for x in lcs if not x.upper().find('NOS')>=0]
    aar=df[lcs].values
    #gmeans=np.exp(np.mean(np.log(aar), axis=1))
    gmeans=np.mean(aar, axis=1)
    gmeans=gmeans.reshape(-1,1)
    #rgmaar=aar/gmeans
    varaar=np.var(aar, axis=1) #change between variance (var) and deviation (std) here
    #mvar=np.mean(varaar, axis=1)
    #mvar=mvar.reshape(-1,1)
    #varaar=varaar/mvar
    varaar=varaar.reshape(-1,1)
    gmaar=np.append(gmeans, varaar, axis=1)
    #gmaar=gmaar[gmaar[:,0].argsort()]
    #sns.set()
    fig=plt.figure(figsize=(12,7))    
    plt.scatter(gmaar[:,0], gmaar[:,1], alpha=0.4)
    sns.regplot(gmaar[:,0], gmaar[:,1], lowess=True, scatter=False, color='r')   
    plt.ylim(gmaar[:,1].min(), gmaar[:,1].max())
    plt.xlim(gmaar[:,0].min(), gmaar[:,0].max())
    plt.xlabel('Mean expression', fontsize=15)
    plt.ylabel('Variance',fontsize=15)
    plt.savefig('norm_var.png')
    plt.close()
    df=pd.DataFrame(aar)
    sns.boxplot(df[df>1.5])
    plt.xlabel('Sample', fontsize=12)
    plt.ylabel('Expresion', fontsize=12)
    #a.set(xlabel='Sample')
    plt.savefig('norm_box.png')
    plt.close()
    print('Figures were generated and saved as norm_var.png. and norm_box.png')
    
def MA_plot(samp1, samp2, name):
    ''' Uses dataframe or array and makes a plot of value difference vs
    mean expression
    name to name a file to save
    eg: df['columnA']
    name is a string variable
    '''
    plt.clf()
    smean=np.mean([samp1.values, samp2.values], axis=0)
    sdif=np.subtract(samp1.values, samp2.values)
    plt.scatter(smean, sdif, alpha=0.4)
    plt.ylim(-2,5)
    plt.xlim(smean.min())
    plt.xlabel('Mean expression')
    plt.ylabel('Difference' )
    plt.tight_layout()
    name='MA_plot'+name+'.png'
    plt.savefig(name)
    plt.close()
    
def expression_plot(samp1, samp2, name):
    '''Uses dataframe or array and makes plot of expression values for two samples
        vs mean expression
    '''
    smean=np.mean([samp1.values, samp2.values], axis=0)
    frames=[x.reshape(-1,1) for x in [samp1, samp2, smean]]
    arr=np.concatenate(frames, axis=1)
    arr=arr[arr[:,2].argsort()]
    print(samp1.name, samp2.name)
    plt.plot(arr[:,0], alpha=0.4, markersize=3, marker='o', 
             color='b', linestyle='none', label=samp1.name)
    plt.plot(arr[:,1], alpha=0.4, markersize=3, marker='o', 
             color='r', linestyle='none', label=samp2.name)
    plt.legend(loc='upper left')
    #plt.ylim()
    ##plt.xlim(smean.min())
    plt.xlabel('Genes')
    plt.ylabel('Expression' )
    name='Expression_plot'+name+'.png'
    plt.savefig(name)
    plt.close()
    
    
def hist_show(sample, df):
#use as a storage for pandas histogram mechanics
    fi=df[sample][df[sample>0.0]].hist(bins=100, range=(0, 1000), alpha=0.4)
    #fi=df[['B18201','E20201']][df['E20201']>0][df['B18201']>0].plot.hist(bins=100, range=(0, 1000), alpha=0.4)
    #pay attention to conditions, do not work as a list, use separately
    #!!! df.hist(bins=40, log=True)
    fi.set_yscale('log')
    plt.savefig('filename.png')
    
    #will continue to reuse the same graph for the next figure untill it closed
    plt.clf()
    # plt.close(fi)

def heat_map(df): #requires #import seaborn as sns#
    plt.clf()
    plt.figure(figsize=(15,15))
    corr=df.corr('pearson')
    sns.heatmap(corr)
    plt.savefig('heat_map.png')
    pass

def wisker():
    pass

def scat_matrix(df, lcs, name='namefig'):
    ''' requres dataframe and list of samples for a scatter. Name for the file to save '''

    plt.clf()
    #does not like NaN values; need to replace and fill before using
    #sns.PairGrid has option dropna. consider
    dfx=df.fillna(0)
    arr=np.array(dfx[lcs])
    arr[arr<0]=0
    dfs=pd.DataFrame(arr, index=dfx.index, columns=lcs)
    #dfx[lcs]=dfx[lcs][dfx[lcs]<0]=0
    sns.set(font_scale=2)
    g=sns.PairGrid(dfs, size=6)
    g=g.map_offdiag(plt.scatter, alpha=0.5)
    g=g.map_diag(plt.hist, bins=50)
    g.set(xlim=(0,6), ylim=(0,6))
    name=name+'png'
    plt.savefig(name)
    plt.close()
    
def cluster_map(dfx):
    plt.clf()
    corr=dfx.corr()
    cmap=sns.diverging_palette(120, 5, as_cmap=True)
    sns.set(font_scale=1.6)
    fig=sns.clustermap(corr, cmap=cmap, linewidths=.5, figsize=(25,25), vmin=0)
    plt.setp(fig.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
    fig.savefig('clustermap.png')


#############MAIN  ###########################
def main():
    
    if len(sys.argv)<2: #check for agrument presence
        print('Usage: graphing.py pickle_df.pkl base_sample \
        substract_sample [substruct_sample2 ...]')
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
    if not isinstance(df, pd.core.frame.DataFrame):
        print('Pickle file is not Pandas DataFrame')
        sys.exit(2)
    
    try:
        del(df['Undetermined'])
    except:
        pass
    
    df=norm_varr(df) #normalize data by total count and log transform
    dfx=anal_prep(df) #set min expression at log(10)=1
    
    if len(sys.argv)>3: ###sss samples are not yet implemented
        qrs=sys.argv[2] #query sample of interest
        sss=sys.argv[3:] #list of samples to substruct
        #call make difference df 
    
        print('Gather 2x fold increased in', qrs, 'over', sss)
        try:
            dff=pd.DataFrame(dfx['Annotation'])
        #dff[qrs]=1
        except:
            pass
        lcf=[x for x in list(dfx.columns) if x!='Annotation' and x!=qrs]
        for i in lcf:
            ij=qrs+'/'+i
            dff[ij]=dfx[qrs]-dfx[i]
        print('dff dataframe contains expression differences')
    
    
if __name__ =='__main__':
    main()
    
