##calculate mean square displacement(MSD) datafiles

import pandas as pd
from ovito.io import import_file, export_file
from ovito.modifiers import *
from ovito.pipeline import *
import numpy as np
from ovito.data import *
from collections import Counter
import os
import glob
import re
import warnings
warnings.filterwarnings('ignore', message='.*OVITO.*PyPI')
#path='/home/heruhe/Desktop/Ga2o3/cascade/5type/1500ev'
#path= '/home/heruhe/Desktop/Ga2o3/cascade/5type_overlapping/cascade_anneal/anneal_dumpfile'
#path='/Users/ruhe/Desktop/ga2o3/FP/tabGap/Ga/anneal_1000k_2ns'
path='/Users/ruhe/Desktop/ga2o3/FP/tabGap/gamma/anneal1000k_2ns'
path='/Users/ruhe/Desktop/ga2o3/overlapping/overlapping_anneal/gamma/2kev_1000k_100ps'
# Use the glob function to search for files that match the pattern
file_pattern = 'anneals1000_100ps.dump*'
matching_files = glob.glob(f"{path}/{file_pattern}")
#print(matching_files)
# Print the list of matching files
pkal=[]
for file in matching_files:
    particles=pd.DataFrame()
    pipeline = import_file(file)
    parts = file.split("/")
    filename = parts[-1]
# Use regular expression to find the number
    match = re.search(r'\d+(?=\D*$)', filename)

    if match:
        i = int(match.group())
    foldern=path+'/1000_100ps_{}'.format(i)
    # Check if the folder exists
    if not os.path.exists(foldern):
    # If it doesn't exist, create the folder
        os.mkdir(foldern)
    pkal.append(i)
# Calculate per-particle displacements with respect to initial simulation frame:
    mod=CalculateDisplacementsModifier()

    pipeline.modifiers.append(mod)
    frame=pipeline.source.num_frames-1

    data = pipeline.compute(frame)
    fn=foldern+'/particle_Displacement_Magnitude.csv'
    particles['Particle Type']=data.particles['Particle Type']
    particles['Particle Identifier']=data.particles['Particle Identifier']
    particles['Displacement Magnitude']=data.particles['Displacement Magnitude']
    particles.to_csv(fn)
pkal=sorted(pkal)
print(pkal)
