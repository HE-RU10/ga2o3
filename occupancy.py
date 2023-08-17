
import pandas as pd
from ovito.io import import_file, export_file
from ovito.modifiers import *
from ovito.pipeline import *
import numpy as np
from ovito.data import *
from collections import Counter
import os


def ovito_analysis(path):
    
  
 

    os.chdir(path) 
    success_s=[]
    for i in range(180):
        if os.path.isfile('./data.lastframe{}'.format(i)):
            success_s.append(i)
    sites=[]
    for i in success_s:
        particles=pd.DataFrame()
        pipeline = import_file("data.lastframe{}".format(i))
        mod=WignerSeitzAnalysisModifier(per_type_occupancies = True)
        mod.reference=FileSource()
        mod.reference.load("data.relaxed")
        pipeline.modifiers.append(mod)
        #pipeline.modifiers.append(ExpressionSelectionModifier(expression = 'Occupancy<=1'))
        #pipeline.modifiers.append(ExpressionSelectionModifier(expression = 'Occupancy.1 + Occupancy.2+Occupancy.3 + Occupancy.4+Occupancy.5 <= 1'))
        #pipeline.modifiers.append(DeleteSelectedModifier())
        #pipeline.modifiers.append(ClusterAnalysisModifier(cutoff=3.425,sort_by_size=True,compute_com=True))
        data = pipeline.compute()
     

        particles['Particle Type']=data.particles['Particle Type']
        #particles['Occupancy']=data.particles['Occupancy']
        #particles['Cluster']=data.particles['Cluster']
        occupancies=pd.DataFrame(data.particles['Occupancy'],columns=['occupancy.1','occupancy.2','occupancy.3','occupancy.4','occupancy.5'])
        particles=pd.concat([particles,occupancies],axis=1)
        particles['frame']=i
        file="PKA_inf{}.txt".format(i)
        f=open(file)
        pkatype=f.readlines()[0].split()[2]

        particles['frame']=i
        particles['pka_type']=pkatype
        particles.to_csv('occupancy_analysis{}.csv'.format(i)) 
        sites.append(particles)
    
    results=pd.concat(sites)
    return results
 
path1='/home/heruhe/Desktop/Ga2o3/cascade/5type/500ev'
ovito_results=ovito_analysis(path1)
ovito_results.to_csv('occupancy_analysis.csv') 
