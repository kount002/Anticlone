# -*- coding: utf-8 -*-
"""
Created on Sat Jun 18 11:28:00 2016
@author: Evgueni
Number of grafs is limited by available colors. Colors are provided by a list in drw_dots function
add more colors is needed. 

"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pickle
import sys
import os

def arr_norm(ll): #normalize array from using total counts, mult by ave(sum)
    ssums=ll.sum(axis=1) #find sums for ea sample (rows)
    ssums=ssums.reshape(ll.shape[0],1) # make a column arr to multiply wt ea row member
    #ssums=ssums.transpose() #
    amp=ll.sum()/ll.shape[0]
    ll=ll*amp/ssums
    return(ll)


def drw_dots(samplist, ll):
    fig=plt.figure(figsize=(7,7)) # adjust figure size
    plt.xlabel(samplist[0], fontsize=15)
    plt.ylabel(samplist[1], fontsize=15)
    plt.xscale('log')
    plt.yscale('log')
    mval=np.amax(ll)
    mval=mval+10
    plt.axis([0.9, mval, 0.9, mval])
    colors=['bo', 'ro', 'go']
    acolors=['blue', 'red', 'green']
    
    for i in range(len(samplist)-1):
        i=i+1
        plt.plot(ll[0,:], ll[i,:], colors[i-1], markersize=3.4)
        po=mval*0.8-(i*0.2*mval)
        plt.annotate(samplist[i], xy=(po,po), xytext=(mval*0.8, po), color=acolors[i-1], fontsize=15)
        
        
    #plt.plot(ll[0,:], ll[1,:], 'bo')
    figname=samplist[0]+'_fig.png'
    figname=os.path.join(inpp, figname)
    print('Figure file name:', figname)
    fig.savefig(figname)
    
def clon_inter(ll): #counts number of intersecting clones that are non-zero in all samples
    clprod=np.prod(ll, axis=0)
    cltot=len(clprod)
    clnz=np.count_nonzero(clprod)
    print('There are ', clnz, 'intersecting clones', round(100*clnz/cltot, 2), '% of total')
    return()

def prune_dic(masdi, r=1): #remove clones with less or equeal then r reads 
    for i in masdi:
        for x,v in masdi[i].items():
           if v[0]<=r:
              masdi[i].pop(x)

if len(sys.argv)<4: #check for agrument presence
    print('Usage: graphing.py pickle.pkl base_sample \
    compare_sample [compare_sample2 ...]')
    sys.exit(2)
    
try:    #check that pickle exists
    with open(sys.argv[1], 'rb') as handle:
        masdi=pickle.load(handle)
except:
    print('Cannot open pickle file. Run clone_count, verify path.')
    sys.exit(2)    

for i in sys.argv[2:]: #check that specified samples are in pickle
    if i not in masdi:
        print ('Specified sample name is not in the pickle file:', i)
        print (masdi.keys())
        sys.exit(2)



masdi=prune_dic(masdi) 

inpp=os.path.dirname(sys.argv[1])
mlist=set() # create master list of all clones
samplist=[]
for i in sys.argv[2:]:
    samplist.append(i) #this to be used for graphs, not here
    #need a set here then convert to list
    
    mlist.update([x for x in masdi[i]])
    #mlist=mlist+[x for x in masdi[i]]   
print('Total clones to graph:', len(mlist))

#mlist=sorted(mlist)
#shlist=[(a,b) for (a, b,c) in mlist]
#mlist=set(shlist)
#print('Was reduced to', len(mlist))    



ll=[]
for i in sys.argv[2:]: # make array of counts out of dictionary
    lis=[masdi[i][x][0] if x in masdi[i] else 0 for x in mlist] #
    ll.append(lis)
    print(i, 'contains', len(masdi[i].keys()), 'clones.')
ll=np.array(ll)
    

ll=arr_norm(ll)    #call normaliztion funcion arr_norm()
clon_inter(ll)
drw_dots(samplist, ll)

