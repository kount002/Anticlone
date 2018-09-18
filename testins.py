
'''
Created on Tue Feb 21 10:23:12 2017
@author: Evgueni
Output a frequency diagram for the length of the clone inserts, 
which determined as a difference in genomic position.
The clones with introns are filtered out. 
uses clone_count dir in .


'''

########### Param ##############
samplenames=['K10120', 'K10320', 'KNOS220', 'K20120', 'K20420', 'K11020', 'K11120']
#not use anymore infile='/clone_count/20118hta.sam'
count=20000 # number of clones to analyse


########### Import #############
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import re

########### Funcions ############

def mk_dic(infile, count):
    print('Collecting on {0}, for {1} lines'.format(infile, count))
    sizelist=[]
    name=''
    with open(infile, 'r') as infl:
        counter=0
        for line in infl: # loop processes single reads, paired see later
            if counter>count:
                break
            lines=line.rstrip().split('\t')
            if lines[7]=='0': # does zero reads first
                continue      # if pair read has value jump ouf of the cycl
            if len(lines[6])>1: # removes reads that al to diff chr (check if sets for names and annotations are present in the list
                continue
            if lines[-1].endswith('not_unique'): #removes lines with multiple alingments
                continue
            if name==lines[0]: #checks if it a reverse read of a pair. Does all math here             
                insize=abs(int(lines[8]))
                cigar1=re.findall(r'M(\d+)N', cigar)
                cigar2=re.findall(r'M(\d+)N', lines[5])
                cigar=list(set(cigar1+cigar2))
                sumcigar=sum([int(x) for x in cigar])
                insertlen=insize-sumcigar
                if insertlen>0: #cannot figure wy some are negative
                    sizelist=sizelist+[insertlen]
                    counter+=1
            name=lines[0]
            cigar=lines[5] 
    return(sizelist)
    


########### Main ################


for i in samplenames:
    samplename=i
    infile='clone_count/'+samplename+'hta.sam'
    sizelist=mk_dic(infile, count)
    sizefilter=[x for x in sizelist if x<700]
    meansize=sum(sizefilter)/len(sizefilter)
    meansize=int(meansize)

    outfile='freqhis_'+samplename+'_'+str(meansize)+'.png'
    print('Mean fragment size:', meansize)
    plt.hist(sizefilter, bins=70)
    plt.ylim(0,4000)
    plt.xlabel('Insert size')
    plt.ylabel('Frequency')
    plt.savefig(outfile, dpi=600)
    plt.close()

