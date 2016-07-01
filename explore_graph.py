''' script to test data using pandas dataframe

'''
import sys
import pandas as pd
import numpy as np
import pickle
import matplotlib
matplotlib.use('Agg')  #need this to turn off figure display
import matplotlib.pyplot as plt




def norm_df(df): 
# clone_count.py outputs raw clone counts. Function normalizes count as: count*ave(sum of sums)/specific sums
    lcs=[x for x in lc if x!='Annotation']
    for col in lcs:
        df[col]=df[col].apply(lambda x: x/df[col].sum()*df[lcs].sum().mean())

    return(df)

















########################## MAIN  ###########################

if len(sys.argv)<4: #check for agrument presence
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



### may not need see and delete
#for i in sys.argv[2:]: #check that specified samples are in pickle
#    if i not in masdi:
#        print ('Specified sample name is not in the pickle file:', i)
#        print (masdi.keys())
#        sys.exit(2)
###

#inpp=os.path.dirname(sys.argv[1])
#mlist=set() # create master list of all clones


df=norm_df(df)



