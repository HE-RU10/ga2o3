
##export xyz file for ovito visualization
import pandas as pd
from ovito.io import import_file, export_file
from ovito.modifiers import *
from ovito.pipeline import *
import numpy as np
from ovito.data import *
from collections import Counter
import os

path='/home/heruhe/Desktop/Ga2o3/cascade/5type/1500ev'
 

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

    export_file(pipeline, "Vsites{}.xyz".format(i), "xyz", columns =["Particle Identifier", "Particle Type", "Position.X", "Position.Y", "Position.Z"])

 



