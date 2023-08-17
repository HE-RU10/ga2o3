
import pandas as pd
from ovito.io import import_file, export_file
from ovito.modifiers import *
from ovito.pipeline import *
import numpy as np
from ovito.data import *
from collections import Counter
from collections import OrderedDict
import os

#calculate interstitial sites and do cluster_analysis using Wigner-Seitz defect analysis function in ovito
def WS_analysis(path):
    a=[]  
    Dict=[]
    ws=pd.DataFrame()
    os.chdir(path) 
    success_s=[]
    for i in range(160):
        if os.path.isfile('./data.lastframe{}'.format(i)):
            success_s.append(i)
    ws['frame']=success_s
    for i in success_s:
        pipeline = import_file("data.lastframe{}".format(i))
        #mod=WignerSeitzAnalysisModifier(output_displaced=True,per_type_occupancies = True)
        mod=WignerSeitzAnalysisModifier(per_type_occupancies = True)
        mod.reference=FileSource()
        mod.reference.load("data.relaxed")


        pipeline.modifiers.append(mod)
        #pipeline.modifiers.append(ExpressionSelectionModifier(expression = 'Occupancy>=1'))
        pipeline.modifiers.append(ExpressionSelectionModifier(expression = 'Occupancy.1 + Occupancy.2+Occupancy.3 + Occupancy.4+Occupancy.5 <= 1'))
        #pipeline.modifiers.append(ExpressionSelectionModifier(expression = 'Occupancy.1 + Occupancy.2+Occupancy.3 + Occupancy.4+Occupancy.5 >= 1'))
        pipeline.modifiers.append(DeleteSelectedModifier())
        pipeline.modifiers.append(ClusterAnalysisModifier(cutoff=3.39,sort_by_size=True,compute_com=True))
        #pipeline.modifiers.append(ClusterAnalysisModifier(cutoff=3.09,sort_by_size=True,compute_com=True))
        data = pipeline.compute()
 
 
        a.append(data.attributes['WignerSeitz.interstitial_count'])
  

    
        cluster_table = data.tables['clusters']
        cls_count=cluster_table['Cluster Size']
        count=Counter(cls_count)
        dictcount=dict(count)
        dictcount=OrderedDict(sorted(dictcount.items()))

        Dict.append(dictcount)
 
    cluster=pd.DataFrame(Dict)
    cluster=cluster.fillna(0)
    cluster['frame']=success_s
    cluster=cluster.set_index('frame')
    ws['interstitial']=a
    ws=ws.set_index('frame')
    
 
    pka=[]
    for i in success_s:
        file="PKA_inf{}.txt".format(i)
        f=open(file)
        pkatype=f.readlines()[0].split()[2]
        pka.append(int(pkatype))
    ws['pka_type']=pka   
    
    results=pd.concat([ws,cluster],axis=1)
    return results
path1='/home/heruhe/Desktop/Ga2o3/cascade/5type/2000ev'
ovito_results=WS_analysis(path1)
for i in range(5):
    i=i+1
    ovito_result=ovito_results[ovito_results['pka_type']==i]
    ovito_result.to_csv('interstitials_analysis{}.csv'.format(i)) 
    #ovito_result.to_csv('vacancy_analysis{}.csv'.format(i)) 

ovito_results.to_csv('interstitials_analysis.csv') 
#ovito_results.to_csv('vacancy_analysis.csv') 
