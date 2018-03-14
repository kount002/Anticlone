# -*- coding: utf-8 -*-
"""
Created on Tue Apr 19 09:42:03 2016
@author: Evgueni
Truncated version of clone_count.py
requires preannotated files
Takes bam file as input and counts all unique clones(read pairs) using unique coordinates

As input the program will require -i
USAGE:
clone-count.py -i input_bam_file [-n output_file_name_prefix] [-b mask_for_batch_processing] [-f bin size]
Example python3.3 clone_count.py -i Processed_data/EK\*/tophat/accepted_hits.bam -f 150 -b yes
Need to escape *
Example nohup python clone_count.py -i recount/K*.sam -f 100 -b yes -n renames
check is done for the length of path(4), processes as a batch
Dictionary structure:
mast:
'Sample 1': ('chr', pos1, pos2):[count, {'Annotation'}]
'Sample 2': ('chr', pos1, pos2):[count, {'Annotation'}]

The script has issue with path expantion when the path contains * symbol, windows and python3.3 clone_count start results in a error message related to unknown arguments. 
either shell or argparse to blame.
However the sripts runs fine if started from IPython3 environment via 'run clone_count' command


"""
########### Param #############
rduce=1 #make 0 if want keep names of the reads for each fragment bin, 1 if want to convert them to counts
maq=10 #removes reads with low mapping quality, mostly ambigious. Set to '0' to disable

###############################

import argparse
import os
import glob
from collections import defaultdict
import sqlite3
import pickle
import pandas as pd
import numpy as np


#import sqlite3
#import csv
#import numpy
#import pysam


def converter(f, sampleout): #converts bam file into sam, annotates it and saves them in subdirecotry (see -n option)
    #pathout=inpp+'/'+sampleout+'hta.sam'
    pathout=os.path.join(inpp, sampleout + 'hta.sam')
    #??Does samtools create new dir as part of output path? May need to make dir with sys.join?
    os.system('samtools sort -on {0} - | \
    samtools view -h - | htseq-count -s no -m intersection-nonempty -o {1} - \
    /nfs/gems_sata/references/hg19/Homo_sapiens/UCSC/hg19/Annotation/Archives/archive-2014-06-02-13-47-56/Genes/genes.gtf \
    >&2'.format(f, pathout))
    #"samtools sort -on Processed_data/"$si"/tophat/accepted_hits.bam - | samtools view -h - "
    return(pathout)
    
def janitor(): #cleanup intermediate files
    pass
    #os.remove() # hta.sam from converter?
    

def collector(pathin, r=10): #process file and gets unique reads
    
    def cust_round(val, r=10): #rounds up the position of the clones to create sorting bins
        rval=int(r*(round(int(val)/r)))
        return(rval)
    
    
    print('Collecting on ', pathin)
    mastdict=defaultdict(list)
    with open(pathin, 'r') as infl:
        for line in infl: # loop processes single reads, paired see later
            lines=line.rstrip().split('\t')
            if lines[7]!='0': # does zero reads first
                continue      # if pair read has value jump ouf of the cycle
            if len(lines[6])>1: # removes reads that al to diff chr (check if sets for names and annotations are present in the list
                continue
            name=[lines[0]]            
            chrn=lines[2]
            pos1=cust_round(lines[3], r)
            pos2=cust_round(lines[7], r)
            #pos1=round(int(lines[3]), r) #convert into int roundup to reduce bins and reduce resolution
            #pos2=round(int(lines[7]), r)  #old function using order 
            annt=[lines[-1].strip('XF:Z:')]
            key=(chrn, pos1, pos2)
            #print(key,'masterkey')
            if len(mastdict[key])==0:   #data structure: {(pos1,pos2):[{set of names},{set of annotations}]}
                mastdict[key].append(set())  #add sets names and annotations to list
                mastdict[key].append(set())
            #print('nameset', name)
            mastdict[key][0].update(name)  #update set of names
            mastdict[key][1].update(annt)
       
    with open(pathin, 'r') as infl:      
        for line in infl: #loop processes paired reads
            lines=line.rstrip().split('\t')
            if lines[7]=='0':
                continue
            if len(lines[6])>1: # removes reads that al to diff chr (check if sets for names and annotations are present in the list
                continue
            if int(lines[4])<maq: # removes reads with mapping quality less then 0.1
                continue
            name=[lines[0]]            
            chrn=lines[2]

            #locate start and end of the clone and orient pos1<pos2. End is start of pos1 plus length of the clone
            pos1=min(int(lines[3]), int(lines[7])) #accouts for orientation of the reads
            pos2=pos1+abs(int(lines[8]))

            #incorrect postion 2 math 
            #pos2=int(lines[7])+int(lines[8])
            #pos1=int(lines[3])
            #pos1=min(pos1, pos2)
            #pos2=max(pos1, pos2)
            
            pos1=cust_round(pos1, r)
            pos2=cust_round(pos2, r)

            ##old solution remove
            ##pos1=cust_round(lines[3], r)
            ##pos2=cust_round(lines[7], r)                        
            ##

            #pos1=round(int(lines[3]), r) #convert into int roundup to reduce bins and reduce resolution
            #pos2=round(int(lines[7]), r) #old 10/100/1000 based rounding
            annt=[lines[-1].strip('XF:Z:')]
            if (chrn, pos2, pos1) in mastdict: #check if mate is present in the mastdict
                continue
            key=(chrn, pos1, pos2)
             
            srv=None
            if (chrn, pos1, 0) in mastdict: #checks and removes single reads
                #print((chrn, pos1, 0), 'zkey')
                srv=mastdict.pop((chrn, pos1, 0)) #saves value of removed SR for future use
                name=name+list(srv[0])
                annt=annt+list(srv[1])
            
            if len(mastdict[key])==0:
                mastdict[key].append(set())
                mastdict[key].append(set())
            mastdict[key][0].update(name)
            mastdict[key][1].update(annt)
             
    return(mastdict)

def mast_mk(pathin): #makes a dict that combines all sample dictionaries 
    mastdict=collector(pathin, r=int(args.bin))
    #saving file to disk
    #? do you need dic.txt if a sql db is created?
    #pathout=inpp+'/'+sample+'dic.txt' #location for the dictionary out file
    #with open(pathout, 'w') as outfl:
    #for k,l in mastdict.items():
    #    print(k,l, file=outfl)
    # '''save file to the direcotry'''
    mast[sample]=mastdict
    return(mast)

def save_db(mast): #(master dict with all clones, ) converts dict with clones into sql db 
    #under development implement clone matching between samples for the same clone    
    ''' one table: (Position-chromosome, position-start, position-end,
    Annotation, Sample1(count), Sample2(count))
    '''
    
    def check(pathdb): #remove db file if present
        if os.path.exists(pathdb):
            conn=sqlite3.connect(pathdb)
            conn.close()
            os.remove(pathdb)
        return()
    #make name for db file
    dbname=inpp+'_hta.db3'
    check(dbname)
    conn=sqlite3.connect(dbname)
    c=conn.cursor()
    
    tbname=inpp+'_tab'
    c.execute('CREATE TABLE {0} ("pos_ch" TEXT, "pos_st" INTEGER, "pos_end" INTEGER, \
    "annotation" TEXT)'.format(tbname))
    #add columns for ea sample
    for sample in mast:
        print('dgs make columns', sample)
        c.execute('ALTER TABLE {0} ADD COLUMN {1} INTEGER'.format(tbname, sample))

    #populate table
    for sample in mast:
        for clone in mast[sample]:
            psch, psst, psend=clone
            annt=', '.join(mast[sample][clone][1])
            count=mast[sample][clone][0]
            
                        
            #print('dgs values', psst, psend, annt, count)
            c.execute('INSERT INTO {0} ("pos_ch", "pos_st", "pos_end", \
            "annotation", {1}) VALUES (?,?,?,?,?)'.format(tbname, sample), \
            (psch, psst, psend, annt, count))
    conn.commit()
    conn.close()
            
def save_pickle(mast): #saves dictionary as a pickle file for use in other scripts
    suff= ('_dic' if type(mast)==dict else '_df')    
    pathf=os.path.join(inpp, inpp+suff+'.pkl')
    with open(pathf, 'wb') as handle:
        pickle.dump(mast, handle)
    return()            

def dic_reduce(mast): #replaces set of read names with count
    for sample in mast:
        for clone in mast[sample]:
            count=len(mast[sample][clone][0])
            mast[sample][clone][0]=count
    return(mast)
    

def dic_df(mast): #convert master dictionary to list of dataframes

    dfmast=pd.DataFrame(index=None, columns=['Annotation'])
    dfmast['Annotation'].astype(str)
    for i in mast:
        dftm=pd.DataFrame.from_dict(mast[i], orient='index')
        dftm.columns=[i, 'Annotation']
        dftm['Annotation']=dftm['Annotation'].astype(str).map(lambda x: x.strip("{'}"))
        dfmast=dfmast.join(dftm, how='outer', lsuffix=i)
        
    lc=['Annotation']+['Annotation'+i for i in mast]
    dfmast[lc]=dfmast[lc].fillna('none')
    #print(dfmast[lc].head())
    dfmast['Annotation']=dfmast[lc].apply(lambda x: tuple(x.unique()), axis=1)
    dfmast['Annotation']=dfmast['Annotation'].apply(lambda x: tuple([c for c in x if c!='none']))
    dfmast.drop(lc[1:], axis=1, inplace=True)
    
    return(dfmast)
    
#def dic_invert(mast): # under development (changes dict structure so the clone is top key)
#    cdic=defaultdict([list])
#    
#    for sample in mast:
#        for clone in mast[sample]:
#            if clone in cdic:
#                pass # present, just add count
#            elif clone: #stopped here
#                pass
#            ch, pos1, pos2 = clone
#            
#            if clone:
#                cdic[clone]
#    
#    return(cdic)
    
            
################# Main ############
       

parser = argparse.ArgumentParser(description='''Converts bam/sam files into dataframe and dictionary for analysis''')
parser.add_argument('-i', '--input', help='input sam/bam file to parse', required=True )
parser.add_argument('-n', '--name_prefix', help='add a prefix to the db name', default='clone_count')
parser.add_argument('-b', '--batch', help='mask to batch files "y"')
parser.add_argument('-f', '--bin', help='bin size for clustering for 10 use "-1", for 100 use "-2" #JUST USE 10, 20, 300 etc.#', default=-1)
args = parser.parse_args()


inp1=args.input

#inp1='clone_count/*.sam' #windos hack for the argparse/special charater

sourcefl=inp1.split('/')
inpp=args.name_prefix
mast={}

print('Source files:', inp1, sourcefl)

if not args.batch: # decide if one file to test or maltiple files to process
    print('not batch processing is in development')
    pass
    
else:
    if sourcefl[-1].endswith(".bam"): # decide if bam file to convert or use already sorted sams files
        try: #try to use the specified bam files
        #if len(sourcefl)==4 and sourcefl[0]=='Processed_data': #check if own processing pipe
            read_fl=glob.glob(inp1)
            os.system('mkdir -p {0}'.format(inpp)) #make output directory
            if not read_fl:
                print('Cannot get any files, check path')   
            for f in read_fl: #starts proceessing of individual files
                sourcefl=f.split('/')
                sample=sourcefl[1]
                print('Processing sample...', sample)
                pathin=converter(f, sample)
                mast=mast_mk(pathin)
        except Exception as e:
            print('Error: ', e)
        #else:
            print('What are the samples names? Diagnostic msg')
    elif sourcefl[-1].endswith('.sam'): # makes sure that we have sam file
        if len(sourcefl)==2:
            read_fl=glob.glob(inp1)
            for f in read_fl: #starts proceessing of individual files
                source=os.path.basename(f) #this is window solution
                sample=source.replace('hta.sam', '', 1)
                #sourcefl=f.split('/')
                #sample=sourcefl[1].strip('hta.sam')
                print('Processing sample...', sample)
                #pathin=inpp+'/'+sample+'hta.sam'
                pathin=os.path.join(inpp, sample+'hta.sam')
                mast=mast_mk(pathin)
                                     
    else:
        print('Cannot figure the path', sourcefl) 
    if rduce!=0:
        dic_reduce(mast)
    #save_db(mast) #convert master dict to sqlite db
    save_pickle(mast) #convert master dict to pickle dump
    save_pickle(dic_df(mast))        
    #dfmast=dic_df(mast)
        
        
             
print('Im done')        

        

#if __name__=='__main__':
#    main()


#def mk_db():
#    inname=(args.name_prefix + 'police.db3')
#    conn=sqlite3.connect(inname)
#    c=conn.cursor()
