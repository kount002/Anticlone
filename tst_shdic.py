########## PARAM ###########
pathin='clone_count/clone_count_dic.pkl'
pathout='clone_count/clone_count_dic_mini.pkl'

######### imports ##########
import pickle






########## main ###########

with open (pathin, 'rb') as f:
    mast=pickle.load(f)

minimast={}
for i in mast.keys():
    minimast[i]={}
    subkey=[x for x in mast[i].keys()]
    print('For {0}, lenth is {1}'.format(i, len(subkey)))
    #minimast[i][subkey[:10]]=mast[i][subkey[:10]]
    for j in subkey[:4000]:
        minimast[i][j]=mast[i][j]
print('Done')
