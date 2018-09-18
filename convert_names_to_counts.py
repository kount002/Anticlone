''' takes non-reduced master.pkl after clone_count.py was ran with rduce=0 parm and converts it to master.pkl with counts'''


########### Param
inpkl='master_names.pkl'  #name of file with names
outpkl='master.pkl'  #name of file with counts (output)



########### Import
import pickle



########## Functions


######### Main

#open non-reduced pickle
with open(inpkl, 'rb') as infl:
    indf=pickle.load(infl)
print('Loaded infile...')

print(len(indf), indf['Annotation'].head(50))

#run reduce function from clone_count
for sample in indf:
    for clone in indf[sample]:
        count=len(indf[sample][clone][0])
        indf[sample][clone][0]=count

print('Converted')

#save a new file
with open(outpkl, 'wb') as outfl:
    pickle.dump(outdf, outfl)

print('Done')


