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
    
    #may need to insert NaN back in place of 0. See what are your needs will be
    #log transfom in the last expression    
    #very fast compared to pandas version '''

    lc=list(df.columns)
    lcs=[x for x in lc if x!='Annotation']
    repls=[x for x in range(10)]  #filter out counts that <10. Do it here to add more weight to highly expressed genes
    df.replace(repls, 0, inplace=True) #filter out counts that <10
    df.fillna(0, inplace=True)
    aar=df[lcs].values
    if method=='med':
        print('Using median normalization.')
        aarm=np.ma.masked_array(aar, mask=aar<=0) #mask 0 counts
        smed=np.ma.median(aarm, axis=0) #median based on columns
        amp=smed.mean()
        aar=(aar*amp)/smed
    else:
        print('Using total count (tc) normalization')
        ssums=aar.sum(axis=0) #sum based on column
        amp=ssums.mean()
        aar=(aar*amp)/ssums

    df[lcs]=pd.DataFrame(aar, index=df.index)
    df.replace(0, np.nan, inplace=True) #replaces all zeros for NaN
    df[lcs]=np.log10(df[lcs]) #log-transfrom data
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
    
