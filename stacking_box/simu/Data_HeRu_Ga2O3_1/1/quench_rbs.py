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
import sys

remote = "yes"

import merRelax as _merr
        
def getFredricCellList(path, character = '*.xyz.gz'):
    """
    Get the list of Fredric's cells and sorted according to the number along the z axis.
    
    path: path to the cell folder;

    Parameters
    ----------
    inpath : TYPE
        DESCRIPTION.
    element : TYPE
        DESCRIPTION.

    Returns
    -------
    None.
    character: the common character of the cells
    
    return:
        the list of the cell names sorted according to the number along the z axis.
    """
    
    file_name = []
    file_number = []
    
    path = path + character
    
    for i in glob.glob(path):
        file_name.append(i)
        file_number.append(i.split('-')[2])
        
    file_number, file_name = zip(* sorted(zip(file_number, file_name) ) )
        
    #for i in range(len(file_name)):
    #    print(file_number[i], file_name[i])
        
    return(file_name)
        
def readFredric(inpath, element):
    """
    Read the cells given by Frederic.
    inpath: path to the cell;
    element: element name.
    """
    
    print("Start to read Frederic's cell: %s\n" % inpath)
    
    coord = []
    eleSequence = []    
    
    with gzip.open(inpath, 'rt') as op:
        lineCount = 0
        for line in op:
            lineCount += 1
            
            if lineCount == 1:
                totNum = int(line)
            elif lineCount == 2:
                if 'Frame' in line and 'number' in line and 'boxsize' in line:
                    linesplit = line.split()
                    boxsize_x = float(linesplit[6])
                    boxsize_y = float(linesplit[7])
                    boxsize_z = float(linesplit[8])
                elif 'Lattice' in line and 'Origin' in line:
                    linesplit = line.split('"')
                    linesplit = linesplit[2].split()
                    
                    boxsize_x = float(linesplit[0])
                    boxsize_y = float(linesplit[4])
                    boxsize_z = float(linesplit[8])
                else:
                    raise KeyError("The format of cell has problem!")
            else:
                linesplit = line.split()
                eleSequence.append(element)
                coord.append([])
                coord[-1].append(float(linesplit[1]))
                coord[-1].append(float(linesplit[2]))
                coord[-1].append(float(linesplit[3]))
                
    coord = np.array(coord)
    eleSequence = np.array(eleSequence)
    
    return(totNum, boxsize_x, boxsize_y, boxsize_z, coord, eleSequence)

def mini_x_y(prog, ctype, flag2, file, eleCol, element, xDataCol, sizex, sizey, sizez):
    """
    Get the minimum x and y latice positions of a cell by using the FCC analyzer.
    Arguments: the arguments of fccAnalyzer
    """
    screen = subprocess.run([prog, ctype, flag2, file, '%s'%eleCol, element,'%s'%xDataCol,\
                             '%s'%sizex, '%s'%sizey, '%s'%sizez],capture_output=True, text=True)
    
    screen = screen.stdout        
    lines = screen.split('\n')
    for i in range(len(lines)):
        
        if "Minimum lattice position" in lines[i]:
            min_x = float(lines[i].split()[3])
            min_y = float(lines[i].split()[4])
            
        if "Reallocating" in lines[i]:
            continue
        elif "Allocating" in lines[i]:
            continue
        else:
            print(lines[i])
    
    return(min_x, min_y)

def writeLmpIn(locout, name, x1, x2, y1, y2, z1, z2, time1=1000, time2=1000):
    """
    Write an input file for LAMMPS.
    locout: output path for the input file;
    name: the name of target cell without the .xyz format and the name for the dump file;
    x1, x2, y1, y2, z1, z2: boundaries on three direactions [A];
    time1: simulation time step for the relax under constant temperature;
    time2: simulation time step for the quench;
    """
    
    with open(locout, "w") as wt:
        wt.write("variable \t name string '%s' \n" % name)
        wt.write("\n")
        wt.write("units \t metal \n")
        wt.write("atom_style	 \t atomic \n")
        wt.write("boundary \t p p p \n")
        wt.write("neighbor \t 2.0 bin \n")
        wt.write("neigh_modify \t delay 0 every 1 check yes \n")
        wt.write("\n")
        #wt.write("region	 \t mybox block %.5f %.5f %.5f %.5f %.5f %.5f units box \n" % \
        #         (x1, x2, y1, y2, z1, z2))
        #wt.write("create_box	 \t 1 mybox\n")
        #wt.write("\n")
        wt.write("log \t log.${name} \n")
        wt.write("\n")
        wt.write("read_data ${name}.xyz \n")
        wt.write("\n")
        wt.write("mass 1 69.723\n")
        wt.write("mass 2 69.723\n")
        wt.write("mass 3 15.9994\n")
        wt.write("mass 4 15.9994\n")
        wt.write("mass 5 15.9994\n")
        wt.write("\n")
        wt.write("pair_style	 \t hybrid/overlay eam/fs tabgap \n")
        if remote == "no":
            wt.write("pair_coeff	 \t * * eam/fs /home/jinxinjx/Licument/Fusion_proA/Data_HeRu_Ga2O3_1/Ga-O_20220921.eam.fs Ga Ga O O O \n")
            wt.write("pair_coeff	 \t * * tabgap /home/jinxinjx/Licument/Fusion_proA/Data_HeRu_Ga2O3_1/Ga-O_20220921.tabgap Ga Ga O O O no yes \n")
        elif remote == "yes":
            wt.write("pair_coeff	 \t * * eam/fs /proj/jinxinjx/MD_Potentials/Ga-O_20220921.eam.fs Ga Ga O O O \n")
            wt.write("pair_coeff	 \t * * tabgap /proj/jinxinjx/MD_Potentials/Ga-O_20220921.tabgap Ga Ga O O O no yes \n")
        else:
            raise ValueError("the remote should be yes or no.")
        wt.write("\n")
        wt.write("thermo	 \t 100 \n")
        wt.write("thermo_style \t custom step temp press pe etotal dt time \n")
        wt.write("\n")
        wt.write("velocity \t all create 290 666 mom yes rot yes dist gaussian \n")
        wt.write("fix \t 1 all npt temp 290 290 0.1 iso 0.0 0.0 0.1 \n")
        wt.write("timestep \t 0.001 \n")
        wt.write("run \t %d\n" % time1)
        wt.write("\n")
        wt.write("unfix \t 1\n")
        wt.write("fix \t 1 all npt temp 290 0.001 0.01 iso 0.0 0.0 0.1\n")
        wt.write("timestep \t 0.001\n")
        wt.write("run \t %d\n" % time2)
        wt.write("write_dump	 \t all custom ${name}.xyz element x y z type modify element Ga1 Ga2 O1 O2 O3\n")
        
def readLmpCell(inpath):
    """
    read the LAMMPS dump .xyz cell

    """
    
    coord = []
    eleSequence = []    
    
    with open(inpath, 'r') as op:
        lineCount = 0
        for line in op:
            lineCount += 1
            linesplit = line.split()
            
            if lineCount == 4:
                totNum = int(line)
            elif lineCount == 6:
                boxsize_x = float(linesplit[1]) - float(linesplit[0])
            elif lineCount == 7:
                boxsize_y = float(linesplit[1]) - float(linesplit[0])
            elif lineCount == 8:
                boxsize_z = float(linesplit[1]) - float(linesplit[0])
            elif lineCount > 9:
                eleSequence.append(linesplit[0])
                coord.append([])
                coord[-1].append(float(linesplit[1]))
                coord[-1].append(float(linesplit[2]))
                coord[-1].append(float(linesplit[3]))
            else:
                continue
                
    coord = np.array(coord)
    eleSequence = np.array(eleSequence)
    
    return(totNum, boxsize_x, boxsize_y, boxsize_z, coord, eleSequence)    
    
        
def t2hms(t):
    """
    Convert time t in seconds to hour-min-s
    """
    time_h = int ( t / 3600 )
    time_min = int (t / 60) - time_h * 60
    time_s = (t) % 60
    
    return(time_h, time_min, time_s)

#=============================================================================#

def quench(path, path2Base, lmp_time1, lmp_time2, nprocessor, lmp_prog):
    """
    The main function for connecting cells.
    Perform a relax after each connection
    
    path: path to cell being relaxed
    path2Base: path to the folder of Base
    deltaZ: the margin distance between two cells
    lmp_time1: the relaxation time step used in LAMMPS simulations
    lmp_time2: the quench time step used in LAMMPS simulations
    nprocessor: number of processors;
    lmp_prog: path to the LAMMPS program

    """
            
    #1. Create cell at1
    at1 = _merr.CELL()
    at1.element = {"Ga1" : 1, "Ga2" : 2, "O1" : 3, "O2" : 4, "O3" : 5}
    at1.mass = {"Ga1" : 69.723, "Ga2" : 69.723, "O1" : 15.9994, "O2" : 15.9994, "O3" : 15.9994}
    at1.massSeqeunce = ["Ga1", "Ga2", "O1", "O2", "O3"]
    at1.readSelf(path)
            
    #2. Write the LAMMPS input file for quench
    lmp_in = "quench.in"
    lmp_out = "tmp"
    writeLmpIn(lmp_in, lmp_out, -1*at1.boxsize_x/2, at1.boxsize_x/2, \
                -1*at1.boxsize_y/2, at1.boxsize_y/2, -1*at1.boxsize_z/2, at1.boxsize_z/2, lmp_time1, lmp_time2)
        
    at1.write_lmp('tmp.xyz')
        
    #4. Relax the cell by using LAMMPS
    subprocess.run(['mpirun', '--mca', 'btl_openib_allow_ib', \
                    '1', '-np', '%s'%nprocessor, lmp_prog, '-in', './%s'%lmp_in])#,capture_output=True, text=True)
        
    #5. Update the at1 cell by reading LAMMPS's dump file
        #1). Update the total number, coordinates and box sizes
    at1.TotNum, at1.boxsize_x, at1.boxsize_y, at1.boxsize_z, at1.xyz, at1.elementSequence = \
        readLmpCell('tmp.xyz')
        #2). Move origin to the center (the origin of the LAMMPS output is very close to the center)
    at1.periodicBC(-at1.boxsize_x/2, at1.boxsize_x/2, -at1.boxsize_y/2, at1.boxsize_y/2, -at1.boxsize_z/2, at1.boxsize_z/2)
        #3). Write the cell
    at1.write('tmp.xyz')
    

#=============================================================================#

time_A = time.time()

print("Runing quench_rbs.py")

"""Initilization parameters"""

#!!!Please pay attention to the path under "remote" (some paths are in the functions)

#Parameters for the quenching
path = ['./final.xyz']
path2Base = '../Base/'
lmp_time1 = 10000
lmp_time2 = 10000
nprocessor = '50'
if remote == 'no':
    lmp_prog = '/home/jinxinjx/Licument/MD_Simu/lammps-2Aug2023/src/lmp_mpi'
elif remote == 'yes':
    lmp_prog = path2Base + 'lmp_mpi'


#Find the target
locin = path

for i in range(len(locin)):  
          
  #Create the target
  quench(path[i], path2Base, lmp_time1, lmp_time2, nprocessor, lmp_prog)
      
  print('Calculating for:', locin[i])


time_end = time.time()

print("Time spent (total): %.3f s"%(time_end - time_A))






