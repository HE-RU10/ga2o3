
##obtain ws defects sites original position vector and type information

import pandas as pd
from ovito.io import import_file, export_file
from ovito.modifiers import *
from ovito.pipeline import *
import numpy as np
from ovito.data import *
from collections import Counter
import os

path='/home/heruhe/Desktop/Ga2o3/cascade/5type/750ev'
 

os.chdir(path) 
success_s=[]
for i in range(180):
    if os.path.isfile('./data.lastframe{}'.format(i)):
        success_s.append(i)
for i in success_s:
    particles=pd.DataFrame()
    pipeline = import_file("data.lastframe{}".format(i))
    mod=WignerSeitzAnalysisModifier(per_type_occupancies = True)
    mod.reference=FileSource()
    mod.reference.load("data.relaxed")
    pipeline.modifiers.append(mod)
    pipeline.modifiers.append(ExpressionSelectionModifier(expression = 'Occupancy.1 + Occupancy.2+Occupancy.3 + Occupancy.4+Occupancy.5 >= 1'))
    pipeline.modifiers.append(DeleteSelectedModifier())
    data = pipeline.compute()
     

    particles['Particle Type']=data.particles['Particle Type']
        #particles['Occupancy']=data.particles['Occupancy']
        #particles['Cluster']=data.particles['Cluster']
    positions=pd.DataFrame(data.particles['Position'],columns=['x','y','z'])
    particles=pd.concat([particles,positions],axis=1)
    file="PKA_inf{}.txt".format(i)
    f=open(file)
    pkatype=f.readlines()[0].split()[2]
    particles['Particle Identifier']=data.particles['Particle Identifier']
    particles['frame']=i
    particles['pka_type']=pkatype
    particles.to_csv('Idefects_original_position{}.csv'.format(i)) 
    particles.to_csv('Vacancy_position{}.csv'.format(i)) 

 



