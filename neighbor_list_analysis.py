import pandas as pd
from ovito.io import import_file, export_file
from ovito.modifiers import *
from ovito.pipeline import *
import numpy as np
from ovito.data import *
from collections import Counter
import os
from ovito.data import CutoffNeighborFinder

path1='/home/heruhe/Desktop/Ga2o3/cascade/5type_overlapping/500ev/out'
#path1='/home/heruhe/Desktop/Ga2o3/cascade/5type/1500ev/anneal'
natom=81920
#natom=160000
os.chdir(path1)
for i in range(1):
    pipeline = import_file("data.lastframe-{}".format(i))
#modifier = CreateBondsModifier(cutoff = 2.3)
#pipeline.modifiers.append(modifier)
    data = pipeline.compute()

# Initialize neighbor finder object:
#cutoff = 2.6


#finder = CutoffNeighborFinder(cutoff, data)
# Set up a neighbor finder for visiting the 12 closest neighbors of each particle.
    finder = NearestNeighborFinder(6, data)

# Prefetch the property array containing the particle type information:
    ptypes = data.particles.particle_types
# Prefetch the property array containing the particle position information:
    position= data.particles.positions
#particle neighbor list
    neighbor=np.zeros((natom,17))
# Loop over all particles:
    for index in range(data.particles.count):
        id=index+1
        neighbor[index,0]=id
        neighbor[index,1] = ptypes[index]
        neighbor[index,2]=position[index][0]
        neighbor[index,3]=position[index][1]
        neighbor[index,4]=position[index][2]

    # Iterate over the neighbors of the current particle:
        j=5
        for neigh in finder.find(index):
        #nid=neigh.index+1
        #print(nid)

        # The index can be used to access properties of the current neighbor, e.g.
            type_of_neighbor = ptypes[neigh.index]
            neighbor[index,j]=type_of_neighbor
            neighbor[index,j+1]=neigh.distance
            j+=2
    ndf=pd.DataFrame(neighbor,columns=['Particle Identifier','Particle Type','x','y','z','n1','d1','n2','d2','n3','d3','n4','d4','n5','d5','n6','d6'])
    ndf.to_csv('ndf_{}.csv'.format(i)) 
