##calculate mean square displacement(MSD) datafiles

import pandas as pd
from ovito.io import import_file, export_file
from ovito.modifiers import *
from ovito.pipeline import *
import numpy as np
from ovito.data import *
from collections import Counter
import os

#path='/home/heruhe/Desktop/Ga2o3/cascade/5type/1500ev'
path= '/home/heruhe/Desktop/Ga2o3/cascade/5type_overlapping/cascade_anneal/anneal_dumpfile'

os.chdir(path) 
particles=pd.DataFrame()

pipeline = import_file("anneals1000_100ps.dump1900")
print(pipeline.source.num_frames)
# Calculate per-particle displacements with respect to initial simulation frame:
mod=CalculateDisplacementsModifier()

pipeline.modifiers.append(mod)
for frame in range(pipeline.source.num_frames):
    data = pipeline.compute(frame)

    particles['Particle Type']=data.particles['Particle Type']
    particles['Particle Identifier']=data.particles['Particle Identifier']
    particles['Displacement Magnitude']=data.particles['Displacement Magnitude']
    particles.to_csv('particle_Displacement_Magnitude{}.csv'.format(frame)) 
