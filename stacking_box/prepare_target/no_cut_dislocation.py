#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 28 10:18:13 2022

@author: jinxinjx

If a dislocation is cut by the z boundaries,
adjust the z boundaries.
"""

import numpy as np
import glob
import subprocess
import os.path

import do_not_cut_dislocation as _noCut

#Find the file name with no 'nocut' inside
file_name = []
for i in glob.glob('./Rotate/data.*'):
    if 'nocut' not in i:
        file_name.append(i)
    
print(file_name)
#Adjust the z boundaries if it is necessary;
#zip the new produced files;
#delete the old files (if a new file is produced)
bin_width = 1
no_cover_limit = 3.1652
for i in range(len(file_name)):
    locin = file_name[i]
    locout = file_name[i]
    
    #_noCut.uncut_Ga2O3(locin, locout, bin_width, no_cover_limit, outType='lammps/data')
    element_Ga2O3 = {1:'Ga1', 2:'Ga2', 3:'O1', 4:'O2', 5:'O3'}
    _noCut.uncut_Ga2O3(locin, locout, bin_width, no_cover_limit, outType='xyz', element=element_Ga2O3)
    
    # if os.path.isfile(locout):
    #     subprocess.run(['gzip', locout])
        
    #     subprocess.run(['rm', locin])
    
    

