import numpy as np
import os
import pandas as pd
path='/home/heruhe/Desktop/Ga2o3/combine_cells/cascade_beta_2box/2000PKA_beta/010'
os.chdir(path) 
atom_number1=81920
box1=pd.read_csv('2000PKA', skiprows=17,sep=' ',nrows=atom_number1,usecols=[1,2,3,4],header=None,names=['type', 'x', 'y', 'z'])
atom_number2=81920
box2=pd.read_csv('ref.dat', skiprows=17,sep=' ',nrows=atom_number2,usecols=[1,2,3,4],header=None,names=['type', 'x', 'y', 'z'])
# Concatenate the two DataFrames vertically
combined_box = pd.concat([box1, box2], ignore_index=True)
# Create a new column 'id' as the index starting from 1
combined_box['id'] = range(1, len(combined_box) + 1)
# Rearrange the columns
combined_box = combined_box[['id', 'type', 'x', 'y', 'z']]
# Specify the file path of the existing file
existing_file_path = 'combined'

# Check if the file exists
try:
    with open(existing_file_path, 'r') as file:
        file.readlines()
except FileNotFoundError:
    # If the file doesn't exist, create it with the header
    combined_box.to_csv(existing_file_path, index=False)
else:
    # If the file exists, append the DataFrame without the header
    combined_box.to_csv(existing_file_path, mode='a', header=False, index=False, sep=' ')





