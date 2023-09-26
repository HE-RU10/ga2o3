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
path= '/home/heruhe/Desktop/Ga2o3/Frankpair/FP_beta/npt/mix/out'

os.chdir(path) 
particles=pd.DataFrame()
success_s=[]
for i in range(1,1001):
    #if os.path.isfile('./data.lastframe{}'.format(i)):
    #    success_s.append(i)
#for s in success_s:
    pipeline = import_file("data.lastframe_{}".format(i))

    # Calculate per-particle displacements with respect to initial simulation frame:
    mod=CalculateDisplacementsModifier()
    mod.reference=FileSource()
    mod.reference.load("data.lastframe_0")
    pipeline.modifiers.append(mod)
    data = pipeline.compute()
    particles['Particle Type']=data.particles['Particle Type']
    particles['Particle Identifier']=data.particles['Particle Identifier']
    particles['Displacement Magnitude']=data.particles['Displacement Magnitude']
    particles.to_csv('particle_Displacement_Magnitude{}.csv'.format(i)) 