''' script to test data using pandas dataframe

'''
import sys
import pandas as pd
import numpy as np
import pickle
import matplotlib
matplotlib.use('Agg')  #need this to turn off figure display
import matplotlib.pyplot as plt


def norm_df(df, lc): # do not use
# clone_count.py outputs raw clone counts. Function normalizes count as: count*ave(sum of sums)/specific sums
    lcs=[x for x in lc if x!='Annotation']
    for col in lcs:
        print('Processing...', col )
        df[col]=df[col].apply(lambda x: x/df[col].sum()*df[lcs].sum().mean())
    return(df)
    
def norm_varr(df):
# clone_count.py outputs raw clone counts. Function normalizes count as: count*ave(sum of sums)/specific sums
#may need to insert NaN back in place of 0. See what are your needs will be
#log transfom in the last expression    
    #very fast compared to pandas version

    lc=list(df.columns)
    lcs=[x for x in lc if x!='Annotation']
    df.fillna(0, inplace=True)
    aar=df[lcs].values
    ssums=aar.sum(axis=0) #sum based on column
    amp=ssums.mean()
    aar=(aar*amp)/ssums
    df[lcs]=pd.DataFrame(aar, index=df.index)
    df.replace(0, np.nan, inplace=True) #replaces all zeros for NaN
    df[lcs]=np.log10(df[lcs]) #log-transfrom data
    
    return(df)

def hist_show(sample, df):
#use as a storage for pandas histogram mechanics
    fi=df[sample][df[sample>0.0]].hist(bins=100, range=(0, 1000), alpha=0.4)
    #fi=df[['B18201','E20201']][df['E20201']>0][df['B18201']>0].plot.hist(bins=100, range=(0, 1000), alpha=0.4)
    #pay attention to conditions, do not work as a list, use separately
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
    dfx=df.fillna(0l)
    arr=np.array(dfx[lcs])
    arr[arr<0]=0
    dfx[lcs]=pd.DataFrame(arr, index=dfx.index)
    #dfx[lcs]=dfx[lcs][dfx[lcs]<0]=0
    sns.set(font_scale=2)
    g=sns.PairGrid(dfx)
    g=g.map(plt.scatter, alpha=0.5)
    plt.savefig('namefig.png')


#############MAIN  ###########################

if len(sys.argv)<2: #check for agrument presence
    print('Usage: graphing.py pickle_df.pkl base_sample \
    compare_sample [compare_sample2 ...]')
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
dfar=norm_varr(df)
