import pandas as pd
from ovito.io import import_file, export_file
from ovito.modifiers import *
from ovito.pipeline import *
import numpy as np
from ovito.data import *
from collections import Counter
from collections import OrderedDict
import os

#calculate interstitial sites recover percentage, eg:data.lastframe9_3_10
def WS_recover(id,shell):
    
    fn='data.lastframe'+str(id)+'_'+str(shell)+'_'
    n=0#total number of data frame
    n_def=0#number of defected frame which FP did notrecover
    for i in range(30):
        fnt='./'+fn+str(i)
        if os.path.isfile(fnt):
            n+=1
            pipeline = import_file(fnt)
            #mod=WignerSeitzAnalysisModifier(output_displaced=True,per_type_occupancies = True)
            mod=WignerSeitzAnalysisModifier(per_type_occupancies = True)
            mod.reference=FileSource()
            mod.reference.load("data.beta_last")


            pipeline.modifiers.append(mod)
            
            data = pipeline.compute()
            if data.attributes['WignerSeitz.interstitial_count']>0:
                n_def+=1 
    return 1-n_def/n
path='/home/heruhe/Desktop/Ga2o3/Frenkelpair/FP_beta/npt/FP_recover/random/O'
#path='/home/heruhe/Desktop/Ga2o3/Frankpair/FP_beta/npt/FP_recover/mix'
os.chdir(path) 
frames=[]
for f in range(1,51):
    fn='data.FP'+str(f)
    pipeline = import_file(fn)
    mod=WignerSeitzAnalysisModifier(per_type_occupancies = True)
    mod.reference=FileSource()
    mod.reference.load("data.beta_last")
    pipeline.modifiers.append(mod)
    data2 = pipeline.compute()
    interstitial_n=data2.attributes['WignerSeitz.interstitial_count']
    if interstitial_n>0:
        frames.append(f)
o_recover=np.zeros((len(frames),4))
for i in range(len(frames)):
    id=frames[i]
    o_recover[i,0]=id
    for s in range(3):
        o_recover[i,s+1]=WS_recover(id,s+1)
print(np.mean(o_recover, axis=0))
o_r=pd.DataFrame(o_recover)
o_r.columns = ['frameid', '1shell', '2shell','3sehll']
o_r.to_csv('o_recover_percentage') 

path='/home/heruhe/Desktop/Ga2o3/Frenkelpair/FP_beta/npt/FP_recover/random/Ga/first_second_shell'
os.chdir(path) 
frames=[]
for f in range(21,71):
    fn='data.FP'+str(f)
    pipeline = import_file(fn)
    mod=WignerSeitzAnalysisModifier(per_type_occupancies = True)
    mod.reference=FileSource()
    mod.reference.load("data.beta_last")
    pipeline.modifiers.append(mod)
    data2 = pipeline.compute()
    interstitial_n=data2.attributes['WignerSeitz.interstitial_count']
    if interstitial_n>0:
        frames.append(f)
Ga_recover=np.zeros((len(frames),3))
for i in range(len(frames)):
    id=frames[i]
    Ga_recover[i,0]=id
    for s in range(2):
        Ga_recover[i,s+1]=WS_recover(id,s+1)
print(np.mean(Ga_recover, axis=0))
Ga_r=pd.DataFrame(Ga_recover)
Ga_r.columns = ['frameid', '1shell', '2shell']
Ga_r.to_csv('ga_recover_percentage') 



  

    

