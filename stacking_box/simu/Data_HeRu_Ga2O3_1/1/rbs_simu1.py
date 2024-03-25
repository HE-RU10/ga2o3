#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 27 15:10:40 2022

@author: jinxinjx
"""

import numpy as np
import gzip
import glob
import subprocess
import time
from shutil import copy2
import os
import sys

remote = "no"

import Create as _create
import Input as _input


def createCellProject5(locin1, pris_x, pris_y, coef_z):
    """
    locin1: the damaged cell
    pris_x, pris_y: the x and y box sizes of a pristine cell
                    (the x and y box sizes of the damaged cell will be adjusted
                     according to these two values)
    coef_z: the coefient of expansion of the damaged cell along the z direction
    """
      
    
    #1. Read the damaged cell
    coord = []
    elementSequence = []
    
    with open(locin1, 'r') as op:
        lineCount = 0
        for line in op:
            lineCount += 1
            linesplit = line.split()
            
            if lineCount == 1:
                TotNum = int(line)
            elif lineCount == 2:
                boxsize_x = float(linesplit[1])
                boxsize_y = float(linesplit[5])
                boxsize_z = float(linesplit[9])
            else:
                elementSequence.append(linesplit[4])
                coord.append([])
                coord[-1].append(float(linesplit[1]))
                coord[-1].append(float(linesplit[2]))
                coord[-1].append(float(linesplit[3]))
                
    W = _create.Crystal()
    
    W.element = {"Ga1" : 1, "Ga2" : 2, "O1" : 3, "O2" : 4, "O3" : 5}
    W.elementSequence = np.array(elementSequence)
    W.TotNum = TotNum
    W.xyz = np.array(coord)
    
    #2. Move the damaged cell origin to zero
    x_min = -1 * boxsize_x / 2 
    y_min = -1 * boxsize_y / 2 
    z_min = -1 * boxsize_z / 2 
    
    W.xyz[:, 0] += (0 - x_min)
    W.xyz[:, 1] += (0 - y_min)
    W.xyz[:, 2] += (0 - z_min)
    
    #3. Adjust the damaged cell size according to the reference sizes
    factor_x = pris_x / boxsize_x 
    factor_y = pris_y / boxsize_y
    factor_z = coef_z
        
    W.xyz[:, 0] *= factor_x
    W.xyz[:, 1] *= factor_y
    W.xyz[:, 2] *= factor_z
    
    boxsize_x *= factor_x
    boxsize_y *= factor_y
    boxsize_z *= factor_z
            
    #6. Write the cell
    with open("coords1.in", 'w') as out:
        out.write("%d\n" % TotNum)
        out.write("Lattice=\"%.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f\"\n" % (boxsize_x,0,0,0,boxsize_y,0,0,0,boxsize_z))
        
        for i in range(len(W.xyz)):
            out.write("%s %.5f %.5f %.5f %s\n"%(W.elementSequence[i],W.xyz[i,0],W.xyz[i,1],W.xyz[i,2],W.element[W.elementSequence[i]]))    
    
    return(boxsize_x, boxsize_y, boxsize_z)
            
#=============================================================================#

time_A = time.time()

"""Initilization parameters"""

#!!!Please pay attention to the path under "remote" (some paths are in the functions)

ener = np.array([1600])

#!!!Consider to change the file names here
nprocessor = '2'
code = 'mpid_rec.1.44'
rbsadec1 = '/home/heruhe/Desktop/Ga2o3/code/ga2o3/stacking_box/simu/Base/' + code
elstop1 = '/home/heruhe/Desktop/Ga2o3/code/ga2o3/stacking_box/simu/Base/He_Ga2O3.in'

#X and y box sizes of pristine cell
pris_x = 97.14422
pris_y = 94.34952
#Coefficient of expansion of the cell along the z direction
coef_z = 1.0

#Special case
simuLength = 10 * 2

dens = 5.668 #[g/cm3] calculated according to the pristine cell

#locin = 'tmp.xyz'
locin = 'final.xyz'
boxSizeX,boxSizeY,boxSizeZ = createCellProject5(locin, pris_x, pris_y, coef_z)

for k in range(len(ener)):

    for j in range(1, 2):
      #Create directory for each energy (top)
      if j == 0:
        direcname = '%d_R' % (ener[k])
        direcin = '%s/in' % (direcname)
        direcout = '%s/out' % (direcname)
        Rotation_theta = 6
        Mode_rotate = 2
      elif j == 1:
        direcname = '%d_A' % (ener[k])
        direcin = '%s/in' % (direcname)
        direcout = '%s/out' % (direcname)
        Rotation_theta = 0
        Mode_rotate = 0
      else:
        raise ValueError('j can only be 0 (random) or 1 (aligned)!')
      
      if not os.path.exists(direcname):
        os.makedirs(direcname)
      else:
        raise ValueError('The directory, called %s, that you want to create exists!'%direcname)
      
      if not os.path.exists(direcin):
        os.makedirs(direcin)
      else:
        raise ValueError('The directory, called %s, that you want to create exists!'%direcin)
    
      if not os.path.exists(direcout):
        os.makedirs(direcout)
      else:
        raise ValueError('The directory, called %s, that you want to create exists!'%direcout)
      #Create directory for each energy (end)
      
      #Move in the base files, and generate other files(top)
      elstop2 = '%s/elstop1.in' % direcin
      rbsadec2 = '%s/' % direcname + code
      
      copy2(elstop1, elstop2)
      copy2(rbsadec1, rbsadec2)
      
      InParam = _input.input_dat()
      ##############################################################################
      ##############################################################################
      InParam.flag = "-rbs"
      
      InParam.Ion = 'He'
      InParam.Energy = ener[k] * 1e3 #[eV]
      InParam.Nchannel = ener[k] / 2
      InParam.Detector_slope = 2 #[keV/channel]
      InParam.Detector_intercept = 1 #[keV]
      InParam.Nhistory = 800000
              
              
      InParam.Theta_min = 0.0
      InParam.Theta_max = 0.12 #IBD characteristic angle
      InParam.Fii_min = 0.0
      InParam.Fii_max = 360.0
      InParam.Rotation_theta = Rotation_theta
      InParam.Rotation_fii = 0.0
      InParam.Mode_rotate = Mode_rotate #0: rotate fii by the user input; 1: rotate fii from 0 to 359 with 1 degree as interval.
      InParam.Mode_ibd = 2
      InParam.Mode_tbm = 0
              
              
      InParam.Ntype = 5
      InParam.Density = dens  #[g/cm3] calculated according to the pristine cell
      InParam.Element = ['GaI', 'GaII', 'OI', 'OII', 'OIII']
      #InParam.Mode_vibrate = 1 # 1:Users input the vibration magnitude
      InParam.Displacement = [0.0625, 0.0625, 0.0723, 0.0723, 0.0723] #[A]
      #InParam.DebyeTemp = 738 #[K] 
      InParam.Temperature = 290 #[K]
      InParam.DebyeDisp = 1 #0 for MD cells
              
              
      InParam.Mini_potential = 0 #[eV]
      InParam.Coll_offset = 0 #[AA]
      InParam.Free_distance = 100 #[AA]
      InParam.Threshold_energy = 1e5 #[eV]
      InParam.Mode_FFP = 1 #0: Poisson process; 1: fixed distance; 2: TRIM
      InParam.Mode_Impact = 1 #0: Iradina; 1: TRIM; 2: Original
      InParam.Elstop_mode = 0
              
              
      InParam.Positionmin_x = 0 #[AA]
      InParam.Positionmin_y = 0 #[AA]
      InParam.Positionmin_z = -3 #[AA]
      InParam.Positionmax_x = boxSizeX #[AA]
      InParam.Positionmax_y = boxSizeY #[AA]
      InParam.Positionmax_z = -3 #[AA]
              
      InParam.Boxsize_x = boxSizeX #[AA]
      InParam.Boxsize_y = boxSizeY #[AA]
      InParam.Boxsize_z = boxSizeZ #[AA]
              
      InParam.Depthsize = simuLength #[AA]
      InParam.DepthStop = simuLength #[AA]
              
      InParam.Boundary_x = 1
      InParam.Boundary_y = 1
      InParam.Boundary_z = 0
              
      InParam.Binsize_x = 3 #[AA]
      InParam.Binsize_y = 3 #[AA]
      InParam.Binsize_z = 3 #[AA]
      
      InParam.DBinsize_x = 3  #Used under "-loc" flag
      InParam.DBinsize_y = 3  #Used under "-loc" flag
      InParam.DBinsize_z = 3  #Used under "-loc" flag
              
              
      InParam.SpreadAngle = 60
      InParam.Seed = 2579
              
              
      InParam.Detector_theta = 165
      InParam.Detector_fii = 0
      InParam.Detector_radius = 10 #[mm]
      InParam.Detector_distance = 5 #[cm]
      InParam.Detector_FWHM = 16.5 #[keV]
              
      if remote == 'no':        
          InParam.table_dat = '/home/heruhe/Desktop/Ga2o3/code/ga2o3/stacking_box/simu/Base/table_beta_Ga2O3.dat'
          InParam.scatter_in = '/home/heruhe/Desktop/Ga2o3/code/ga2o3/stacking_box/simu/Base/Scatter.mat'
          InParam.time_in = '/home/heruhe/Desktop/Ga2o3/code/ga2o3/stacking_box/simu/Base/Time.mat'
      else:
          InParam.table_dat = '../../Base/table_beta_Ga2O3.dat'
          InParam.scatter_in = '/proj/jinxinjx/RBSADEC/Data/Scatter.mat'
          InParam.time_in = '/proj/jinxinjx/RBSADEC/Data/Time.mat'
      
      InParam.elstop1_in = './in/elstop1.in'
      InParam.coords1_in = '../coords1.in'
      
  
      InParam.regionBox = [0, 0, 0, 0, 0, 0, 0, 0] #Used under "-loc2" flag
      InParam.regionSphere = [0, 0, 0, 0, 0, 0] #Used under "-loc2" flag
      ##############################################################################
      ##############################################################################
      
      InParam.write_dat("%s/input.dat" % direcin)
      #Move in the base files, and generate other files(end)
      
      
      #Run the RBSADEC code
      runCode = 'mpirun --mca btl_openib_allow_ib 1 -np ' + nprocessor + ' ./' + code
      os.chdir('./%s' % direcname) # Go to the sub directory where the RBSADEC code locates
      cmd='%s %s'%(runCode, InParam.flag)
      os.system(cmd)
      #Come back to the location of this python file
      os.chdir('..') 
      
      if os.path.isfile("%s/elstop1.in" % direcin):
          cmd = 'rm %s/elstop1.in' % direcin
          os.system(cmd)
      if os.path.isfile("%s/" % direcname + code):
          cmd = 'rm %s/' % direcname + code
          os.system(cmd)
      print('=======BRAVO, SIMULATION IN %s IS FINISHED!=======\n' % (direcname))


if os.path.isfile("%s/../../coords1.in" % direcin):
   cmd = 'rm %s/../../coords1.in' % direcin
   os.system(cmd)



time_end = time.time()

print("Time spent (total): %.3f s"%(time_end - time_A))






