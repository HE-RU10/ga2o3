#neighbor_list analysis defects
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os 

path1='/home/heruhe/Desktop/Ga2o3/cascade/5type_overlapping/500ev/out'
natom=81920
os.chdir(path1)
for f in range(1):
    #read 6 nearest neighbor info into ndf dataframe
    ndf=pd.read_csv('ndf_{}.csv'.format(f))
    ndf.drop(ndf.columns[0],axis=1,inplace=True)
    ##Ga1 perfect crystal
    nga1=ndf[ndf['Particle Type']==1] #all Ga1 df
    nga1['nsum4']= nga1.apply(lambda row: row['n1'] + row['n2']+row['n3']+row['n4'], axis=1)
    nga1['nsum6']= nga1.apply(lambda row: row['n1'] + row['n2']+row['n3']+row['n4']+row['n5']+row['n6'], axis=1)
    npga1=nga1[(nga1['nsum4']==16)|(nga1['nsum6']==25)] #perfect Ga1 crystal
    ##Ga2 perfect crystal
    nga2=ndf[ndf['Particle Type']==2] #all Ga2 df
    nga2['nsum4']= nga2.apply(lambda row: row['n1'] + row['n2']+row['n3']+row['n4'], axis=1)
    nga2['nsum6']= nga2.apply(lambda row: row['n1'] + row['n2']+row['n3']+row['n4']+row['n5']+row['n6'], axis=1)
    npga2=nga2[(nga2['nsum4']==16)|(nga2['nsum6']==25)] #Ga2 perfect crystal df

    #Ga1 defects
    ndga1=nga1[~nga1.isin(npga1)].dropna()#Ga1 defects
    dropl=[]
    for i in ndga1['Particle Identifier']:
        x=ndga1[ndga1['Particle Identifier']==i]['x'].values[0]
        z=ndga1[ndga1['Particle Identifier']==i]['z'].values[0]
        n=len(npga2[(abs(npga2['x']-x)<0.9)&(abs(npga2['z']-z)<0.9)])+len(npga1[(abs(npga1['x']-x)<0.9)&(abs(npga1['z']-z)<0.9)])  
        if n>2:    
            dropl.append(int(i-1))
    ndga1 = ndga1.drop(index=dropl)

    #Ga2 defects
    ndga2=nga2[~nga2.isin(npga2)].dropna()#Ga2 defects
    dropl=[]
    for i in ndga2['Particle Identifier']:
        x=ndga2[ndga2['Particle Identifier']==i]['x'].values[0]
        z=ndga2[ndga2['Particle Identifier']==i]['z'].values[0]
        n=len(npga2[(abs(npga2['x']-x)<0.9)&(abs(npga2['z']-z)<0.9)])+len(npga1[(abs(npga1['x']-x)<0.9)&(abs(npga1['z']-z)<0.9)])    
        if n>2:
            dropl.append(int(i-1))
    ndga2 = ndga2.drop(index=dropl)

    
    defects=pd.concat([ndga1,ndga2])
    # prepare output folder
    if not os.path.isdir('frame{}'.format(f)):
        os.mkdir('frame{}'.format(f))
    defects.to_csv('frame{}/defects{}.csv'.format(f,f))     


