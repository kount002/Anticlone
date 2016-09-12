''' script to test data using pandas dataframe

'''
import sys
import pandas as pd
import numpy as np
import pickle
import seaborn as sns
#import matplotlib
#matplotlib.use('Agg')  #need this to turn off figure display
import matplotlib.pyplot as plt


def norm_df(df, lc): # do not use this function too slow
# clone_count.py outputs raw clone counts. Function normalizes count as: count*ave(sum of sums)/specific sums
    lcs=[x for x in lc if x!='Annotation']
    for col in lcs:
        print('Processing...', col )
        df[col]=df[col].apply(lambda x: x/df[col].sum()*df[lcs].sum().mean())
    return(df)
    
def norm_varr(df, method='tc'):
    ''' clone_count.py outputs raw clone counts. 
    Takes values for normalization methods: 'tc' - total count or 'med' - median
    'tc': (default) Function normalizes count as: count*ave(sum of sums)/specific sums
    'med':  Function normalizes count as: count*ave(of medians)/specific median
    does ignores 0 counts for median calculation.
    'max': uses max value to normalize
    'upper': uses upper quartile normalization
    Does not consider samples with 'NOS' in name when removes low count genes from analysis
    
    #may need to insert NaN back in place of 0. See what are your needs will be
    #log transfom in the last expression    
    #very fast compared to pandas version '''

    lc=list(df.columns)
    lcs=[x for x in lc if x!='Annotation']
    lcsnos=[x for x in lcs if not x.upper().find('NOS')>=0] #exclude NOS samples (for reduction)
    #repls=[x for x in range(10)]  #filter out counts that <10. Do it here to add more weight to highly expressed genes
    #df.replace(repls, 10, inplace=True) #filter out counts that <10
    #df=df.loc[(df[lcs]>0).all(axis=1),:]
    df.fillna(0, inplace=True)
    aar=df[lcs].values
    aar[aar<2]=0
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
        print(supper, perc)
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
        meansample=np.exp(np.mean(np.log(aar), axis=1))
        meansample=meansample.reshape(-1,1)
        ratsample=aar/meansample
        supper=np.percentile(ratsample, perc, axis=0)
        amp=np.exp(np.mean(np.log(supper)))
        supper=supper.reshape(1,-1)
        aar=(amp*aar)/supper
    else:
        pass #for none normalization

    df[lcs]=pd.DataFrame(aar, index=df.index)
    df[df<10]= 10 #filter out counts that <10
    df.replace(0, np.nan, inplace=True) #replaces all zeros for NaN #doesn't make sence anymore
    df[lcs]=np.log10(df[lcs]) #log-transfrom data
    #df=df.loc[(df[lcs]>1).any(axis=1),:]
    df=df.loc[(df[lcsnos]>1).any(axis=1),:] #may adjust value for a filter
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

def norm_plot(df):
    ''' plot variance of all samples vs mean expression
        helps to determine which normalization works best for a set 
    '''
    
    lc=list(df.columns)
    lcs=[x for x in lc if x!='Annotation']
    aar=df[lcs].values
    #gmeans=np.exp(np.mean(np.log(aar), axis=1))
    gmeans=np.mean(aar, axis=1)
    gmeans=gmeans.reshape(-1,1)
    #rgmaar=aar/gmeans
    varaar=np.std(aar, axis=1)
    
    #mvar=np.mean(varaar, axis=1)
    #mvar=mvar.reshape(-1,1)
    #varaar=varaar/mvar
    varaar=varaar.reshape(-1,1)
    gmaar=np.append(gmeans, varaar, axis=1)
    #gmaar=gmaar[gmaar[:,0].argsort()]
    
    #sns.set()    
    plt.scatter(gmaar[:,0], gmaar[:,1], alpha=0.4)
    sns.regplot(gmaar[:,0], gmaar[:,1], lowess=True, scatter=False, color='r')
    plt.ylim(0, gmaar[:,1].max())
    plt.xlim(1, gmaar[:,0].max())
    plt.xlabel('Mean expression')
    plt.ylabel('Standard deviation')
    plt.savefig('norm_var.png')
    plt.close()
    df=pd.DataFrame(aar)
    sns.boxplot(df[df>1.5])
    plt.xlabel('Sample', fontsize=12)
    plt.ylabel('Expresion', fontsize=12)
    #a.set(xlabel='Sample')
    plt.savefig('norm_box.png')
    print('Figures were generated and saved as norm_var.png. and norm')
    
def MA_plot(samp1, samp2, name):
    ''' Uses dataframe or array and makes a plot of value difference vs
    mean expression
    name to name a file to save
    '''
    smean=np.mean([samp1.values, samp2.values], axis=0)/2
    sdif=np.subtract(samp1.values, samp2.values)
    plt.scatter(smean, sdif, alpha=0.4)
    #plt.ylim()
    plt.xlim(smean.min())
    plt.xlabel('Mean expression')
    plt.ylabel('Difference' )
    name='MA_plot'+name+'.png'
    plt.savefig(name)
    

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

def heat_map(): #requires #import seaborn as sns#
    plt.clf()
    corr=df[df[lcs]>10].corr('pearson')
    sns.heatmap(corr)
    plt.savefig('heat_map.png')
    pass

def wisker():
    pass

def scat_matrix(): #need lcs
    plt.clf()
    #does not like NaN values; need to replace and fill before using
    #sns.PairGrid has option dropna. consider
    dfx=df.fillna(0)
    arr=np.array(dfx[lcs])
    arr[arr<0]=0
    dfx[lcs]=pd.DataFrame(arr, index=dfx.index)
    #dfx[lcs]=dfx[lcs][dfx[lcs]<0]=0
    sns.set(font_scale=2)
    g=sns.PairGrid(dfx)
    g=g.map(plt.scatter, alpha=0.5)
    plt.savefig('namefig.png')
    
def cluster_map():
    plt.clf()
    corr=dfx.corr()
    cmap=sns.diverging_palette(120, 5, as_cmap=True)
    sns.set(font_scale=1.6)
    fig=sns.clustermap(corr, cmap=cmap, linewidths=.5)
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
    
