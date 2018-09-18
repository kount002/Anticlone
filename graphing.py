# -*- coding: utf-8 -*-
"""
Created on Sat Jun 18 11:28:00 2016
@author: Evgueni
Number of grafs is limited by available colors. Colors are provided by a list in drw_dots function
add more colors is needed. 

"""

import matplotlib.pyplot as plt
import numpy as np
import pickle
import sys

def arr_norm(ll): #normalize array from using total counts, mult by ave(sum)
    ssums=ll.sum(axis=1) #find sums for ea sample (rows)
    ssums=ssums.reshape(ll.shape[0],1) # make a column arr to multiply wt ea row member
    #ssums=ssums.transpose() #
    amp=ll.sum()/ll.shape[0]
    ll=ll*amp/ssums
    return(ll)


def drw_dots(samplist, ll):
    fig=plt.figure(figsize=(7,7)) # adjust figure size
    plt.xlabel(samplist[0], fontsize=14)
    #plt.xscale('log')
    #plt.yscale('log')
    mval=np.amax(ll)
    mval=mval+1
    plt.axis([-1, mval, -1, mval])
    colors=['bo', 'ro', 'go']
    acolors=['blue', 'red', 'green']
    
    for i in range(len(samplist)-1):
        i=i+1
        plt.plot(ll[0,:], ll[i,:], colors[i-1])
        po=60-i*5
        plt.annotate(samplist[i], xy=(2,1), xytext=(60, po), color=acolors[i-1], fontsize=14)
        
        
    #plt.plot(ll[0,:], ll[1,:], 'bo')
    
    fig.savefig('filename.pdf')
    

#plt.plot(s1list, s2list, 'bo', s1list, s3list, 'ro', markersize=4)
#x=[1,2,3,4] 
#y=[20, 25, 27, 29]
#
#t=np.arange(0., 5., 0.2)
#r=np.arange(0., 5., 0.1)
#
#plt.figure(figsize=(9,9))
#plt.subplot(221)
#plt.plot(t, t**2, 'bo', r, r**3, 'k', markersize=.3)
#
#plt.subplot(222)
#plt.plot(t, t**2, 'bo', r, r**3, 'r^', markersize=4)

if len(sys.argv)<3: #check for agrument presence
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
    if not masdi[i]:
        print ('Specified sample name is not in the pickle file:', i)
        print (masdi.keys())
        sys.exit(2)


mlist=[] # create master list of all clones
samplist=[]
for i in sys.argv[2:]:
    samplist.append(i)
    mlist=mlist+[x for x in masdi[i]]   
print('Total clones to graph:', len(mlist))


ll=[]
for i in sys.argv[2:]: # make array of counts out of dictionary 
    lis=[masdi[i][x][0] if x in masdi[i] else 0 for x in mlist]
    ll.append(lis)
ll=np.array(ll)
    
#ll=np.empty([1,len(mlist)], dtype=float) #create list of lists for count in samples. first is base
#for i in sys.argv[2:]:
#    lis=[masdi[i][x][0] if x in masdi[i] else 0 for x in mlist]
#    lis=np.array(lis)
#    np.vstack([ll, lis])
    #ll.append(lis)
    

ll=arr_norm(ll)    #call normaliztion funcion arr_norm()
drw_dots(samplist, ll)



#f='clone_count.pkl'
#s1='A16000'
#s2='B18201'
#s3='C18RUT'
#
#with open(f, 'rb') as handle:
#    mas=pickle.load(handle)
#comb=set()
#comb=set(mas[s1])
#
##print('s1', len(comb))
#comb.update(mas[s2])
#comb.update(mas[s3])
##print('s1 and s2', len(comb))
#
#s1list=[]
#s2list=[]
#s2list=[]
#
#comb=list(comb)
#s1list=[mas[s1][i][0] if i in mas[s1] else 0 for i in comb]
#s2list=[mas[s2][i][0] if i in mas[s2] else 0 for i in comb]
#s3list=[mas[s3][i][0] if i in mas[s3] else 0 for i in comb]
#
#fig=plt.figure(figsize=(7,7))
#plt.plot(s1list, s2list, 'bo', s1list, s3list, 'ro', markersize=4)
##plt.xscale('log')
##plt.yscale('log')
#plt.xlabel(s1)
#
#allval=s1list+s2list+s3list
#amax=max(allval)
#plt.axis([-1, amax+1, -1, amax+1])
#
#fig.savefig('test.pdf')





#fig=plt.plot(t, t, 'r--', t, t**2, 'bs', t, t**3, 'g^', linewidth=1)
#fig=plt.ylabel('yyy')


#plt.axis([0,6,15,30])
#plt.show(fig)