'''
script is designed to identify annotated sequences.
Takes a pkl output of clone_count.py with names that were not reduced to the pure counts (1)
Takes also a csv file that is produced by clone_meandiff.py with not collapsed coordinates that still have start and end.
based on coordinates the script removes short bins less than 150 and makes a ditionary of the coordinates and names


'''


############ param ##########

indf='fold_out_master.csv' #CSV that lists unannotated fragments
indic='test_sample/clone_count/clone_count_dic.pkl' #pickle with dictionatry that contains names and not counts, turn off counting function( dic_reduce) in clone_count.py

fltrshort=200 # set to remove every clone that is shorter 
toplist=100 # set number of frabment to investigate
sambamloc='clone_count/*.sam' #location of sam or bam files


############# imports ########
import pickle
import pandas as pd
import explore_graph as exg
import glob
import pandas as pd


############# functions ######
def filter_short(df):
    def strips(x): #removes simbols from index value
        for i in ['"', "'", '(', ')', ' ']:
            x=x.replace(i, '')
        x=tuple(x.split(','))
        return(x)
    '''remove short fragments and single reads '''
    loc=df.index
    loc=loc.map(strips)
    df['Chr'], df['Pos1'], df['Pos2']=zip(*loc)
    df[['Pos1', 'Pos2']]=df[['Pos1', 'Pos2']].astype('int')
    print('Length of array before filtering short fragment', len(df))
    df=df[abs(df['Pos1']-df['Pos2'])>fltrshort] #dsg: check if int
    print('Length of array after filtering fragments shorter than', fltrshort, len(df))
    df=df[(df['Pos2']!=0) & (df['Pos2']!=0)] #removes bins with unpaired ends
    subloc=df[['Chr', 'Pos1','Pos2']]
    sloc=[tuple(x) for x in subloc.values]
    print('Length of array after filtering unpaired bins', len(df))
    #df.set_index([['Chr','Pos1','Pos2']], inplace=True)
    #print(df.columns)
    del(df['Chr'])
    del(df['Pos1'])
    del(df['Pos2'])
    return(df, sloc)


def blast_seq(seq):
    from Bio.Blast import NCBIWWW
    from Bio.Blast import NCBIXML

    b_result=NCBIWWW.qblast('blastn', 'nt', seq, descriptions=1, hitlist_size=1)
    b_parsed=NCBIXML.parse(b_result)
    item=next(b_parsed)
    for i in item.alignments:
        titl=i.title
        print(titl)
        return(titl) 




def extractor(nm, sambamloc):
    '''extracts a sequence read form bam file (sambamloc) using a read name(nm)'''
    readsam=glob.glob(sambamloc)
    for i in readsam:
        if i.endswith('bam'):
            print('Cannot read bam files yet')
            #open convertion form bam to sam
        #read sam
        print('Looking for sequence in {0} file'.format(i))
        with open(i, 'r') as samfile:
            for line in samfile:
                liner=line.rstrip().split('\t')
                if nm==liner[0]:
                    print('Blasting', liner[0], liner[9])
                    #start blast module, break  out of cycle (test only one read)
                    titl=blast_seq(liner[9])
                    return (titl)



############ main ###########

#create a list of coordinates to extract
df=pd.read_csv(indf, index_col='Unnamed: 0') #open csv with counts 
df, sloc=filter_short(df)   #filter shorts
#df=df.iloc[:toplist, :]
#srch=list(df.index)
'''discard df?'''
sloc=sloc[:toplist]

#open dictionary with names and extract sequences 
mast=exg.open_lib(indic)
anndict={}
for i in sloc: #takes a coordinate
    for j in mast.keys(): #takes a sample
        names=mast[j][i]
        if names:
            print('Found match in master dictionary:', j, i, names[0])
            #extract sequences by names
            namelist=list(names[0])
            #anndict={}
            for n in namelist:
                titl=extractor(n, sambamloc) # run function to extract sequece 
                si=str(i)
                anndict[si]=titl
                with open('out.txt', 'a') as fil:
                    line=si+','+titl+'\n'
                    fil.write(line)
                break
            break

dfann=pd.DataFrame.from_dict(anndict, orient='index')
dfann.columns=['Blast_annotation']
dfall=df.merge(right=dfann, how='outer', left_index=True, right_index=True)
outdf='blast_'+indf
dfall.to_csv(outdf)

print('done')
