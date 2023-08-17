import numpy as np
import subprocess
import ase.io
from ase.data import atomic_numbers
import os
import sys
#transfer lammps data to ase sorted data
path1='/home/heruhe/Desktop/Ga2o3/cascade/5type_overlapping/500ev/out'
#path1='/home/heruhe/Desktop/Ga2o3/cascade/5type/1500ev/anneal'
os.chdir(path1)

datafile = 'data.lastframe0'  # input data file, should be relaxed lammps-data with velocities!

elements = ['B','Be','H','He','Li']  # to keep atom type number,assign type 1 to B,2 to Be, 3 to H, 4 to He, 5 to Li
elements_real = ['O1','O2','O3','Ga2','Ga1'] 
# read in data correctly... map Z and sort by id
Zs = [atomic_numbers[s] for s in elements]
Zs.insert(0, 0)  # 1st lammps type = 1
atoms = ase.io.read(datafile, format='lammps-data', style='atomic', Z_of_type=Zs,sort_by_id=True)
# save output
ase.io.write(f'data.lastframe-0', atoms, format='lammps-data',velocities=True)
    
   