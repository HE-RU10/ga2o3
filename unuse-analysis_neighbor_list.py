#atom's neighbor_list
import pandas as pd
import numpy as np
import os 


#obtain defects atoms' neighbor list, combine the list into defects atoms dataframe
def calculate_neighbor(I_atomsdf,df,d,cutoff=2.3):
    #input:
    #I_atomsdf is defects atoms dataframe,including position & index& type info
    #df is all atoms dataframe,,including position & index& type info
    #d = width/2 of the cube choosing to calculate center atoms neighbor list
    # output:
    # neighbordf is a dataframe ,including  index& type &position & 1-6 neighbor list  info
    neighbors=np.zeros((len(I_atomsdf['Particle Identifier']),17))
    for i in range(len(I_atomsdf['Particle Identifier'])):
        index=I_atomsdf['Particle Identifier'].tolist()[i]
        #obtain interstitial atom's positon
        x=I_atomsdf[I_atomsdf['Particle Identifier']==index]['x'].values
        y=I_atomsdf[I_atomsdf['Particle Identifier']==index]['y'].values
        z=I_atomsdf[I_atomsdf['Particle Identifier']==index]['z'].values

        #regard the interstitial atom as center, choose a cube distrcition, extract all atoms info in the cube  
        cube_df = df[(np.abs(df.x-x)<d)& (np.abs(df.y-y)<d)&(np.abs(df.z-z)<d)]
        #In the cube, calculate each atom's distance with the center interstitial atom
        cube_df.loc[:,'distance']=np.sqrt((cube_df.x-x)**2+(cube_df.y-y)**2+(cube_df.z-z)**2)
        cube_df.sort_values(by=['distance'],inplace=True)
        neighbors[i,0]=I_atomsdf[I_atomsdf['Particle Identifier']==index]['Particle Identifier'].values
        neighbors[i,1]=I_atomsdf[I_atomsdf['Particle Identifier']==index]['Particle Type'].values
        neighbors[i,2]=x
        neighbors[i,3]=y
        neighbors[i,4]=z
        #extract 6 nearest atoms' particle type info to obtain neighbor list,cutoff distance is 2.2
    
        for d in range(1,7):

            if cube_df['distance'].values[d]<=cutoff:
                neighbors[i,d+4]=cube_df['Particle Type'].values[d]
                neighbors[i,d+10]=cube_df['distance'].values[d]
               
            else:
                continue
               
    neighbordf=pd.DataFrame(neighbors)
    neighbordf.columns=['Particle Identifier','Particle Type','x','y','z','n1','n2','n3','n4','n5','n6','d1','d2','d3','d4','d5','d6']
    return  neighbordf
    
    
    
path1='/home/heruhe/Desktop/Ga2o3/cascade/5type_overlapping/out'
natom=81920
os.chdir(path1)
df5=pd.read_csv('data.lastframe-5', skiprows =11,header=None, nrows=natom,sep="\s+")


df5.columns=['Particle Identifier','Particle Type','x','y','z']

ndf5=calculate_neighbor(df5,df5,4)
ndf5.to_csv('ndf5.csv') 