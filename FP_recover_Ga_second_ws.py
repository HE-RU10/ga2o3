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
def WS_recover(id):

# Specify the directory path
    folder_path = "/home/heruhe/Desktop/Ga2o3/Frenkelpair/FP_beta/npt/FP_recover/random/Ga/second_shell/second_shell/"+str(id)

# List all files in the specified directory
    files = os.listdir(folder_path)
    os.chdir(folder_path) 
# Filter files that start with 'data.last'
    matching_files = [file for file in files if file.startswith("data.lastframe")]
    n=0#total number of data frame
    n_def=0#number of defected frame which FP did notrecover
    for fn in matching_files:
            pipeline = import_file(fn)
            #mod=WignerSeitzAnalysisModifier(output_displaced=True,per_type_occupancies = True)
            mod=WignerSeitzAnalysisModifier(per_type_occupancies = True)
            mod.reference=FileSource()
            mod.reference.load("/home/heruhe/Desktop/Ga2o3/Frenkelpair/FP_beta/npt/FP_recover/random/Ga/first_second_shell/data.beta_last")
            pipeline.modifiers.append(mod)
            
            data = pipeline.compute()
            if data.attributes['WignerSeitz.interstitial_count']==0:
                n_def+=1 
    return (n_def, len(matching_files))

path='/home/heruhe/Desktop/Ga2o3/Frenkelpair/FP_beta/npt/FP_recover/random/Ga/second_shell/second_shell'
frames=os.listdir(path)

Ga_recover=np.zeros((len(frames),3))
j=0
for i in frames:
    id=int(i)
    Ga_recover[j,0]=i
    
    Ga_recover[j,1]=WS_recover(id)[0]
    Ga_recover[j,2]=WS_recover(id)[1]
    j+=1
Ga_r=pd.DataFrame(Ga_recover)
Ga_r.columns = ['frameid', 'recover', 'total_amount']
Ga_r.to_csv('/home/heruhe/Desktop/Ga2o3/Frenkelpair/FP_beta/npt/FP_recover/random/Ga/second_shell/second_shell/ga_secondshell_recover_percentage') 



  
