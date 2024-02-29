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

remote = "yes"


class CELL:
    """
    """
    
    def __init__(self):
        
        self.element = {"H": 1}
        self.elementSequence = []
        
        self.xyz = np.empty([1,3])
        self.TotNum = 0
        self.boxsize_x, self.boxsize_y, self.boxsize_z = 2, 2, 2
        
        self.min_x, self.min_y = -1, -1
        
        self.centerFlag = 0
        
    def translate(self, x, y, z):
        """
        Translate the cell along three directions by distances x, y and z
        """
        
        if (x != 0):
            self.xyz[:, 0] += x
        if (y != 0):
            self.xyz[:, 1] += y
        if (z != 0):
            self.xyz[:, 2] += z
            
        self.centerFlag = 0
            
    def shrink(self, x, y, z):
        """
        Shrink / expand the cell along the three direction by factors of x, y and z
        """
        
        if (x != 1):
            self.xyz[:, 0] *= x
            self.boxsize_x *= x
        if (y != 1):
            self.xyz[:, 1] *= y
            self.boxsize_y *= y
        if (z != 1):
            self.xyz[:, 2] *= z
            self.boxsize_z *= z
            
        self.centerFlag = 0
            
    def centerOrigin(self):
        """
        Move the origin from the corner to the cell center
        """
        
        self.xyz[:, 0] -= (self.boxsize_x / 2)
        self.xyz[:, 1] -= (self.boxsize_y / 2)
        self.xyz[:, 2] -= (self.boxsize_z / 2)
        
        self.centerFlag = 1
        
    def periodicBC(self, x1, x2, y1, y2, z1, z2):
        """
        Perform a periodic boundary condtion on cells.
        The boundary limits: x1, x2, y1, y2, z1, z2
        
        (Note, the function can only handle atoms in the first neighbour cells)
        """
        
        limit_up = [x2, y2, z2]
        limit_down = [x1, y1, z1]
        move = [self.boxsize_x, self.boxsize_y, self.boxsize_z]
        
        for i in range(self.TotNum):
            for j in range(3):
                if self.xyz[i, j] >= limit_up[j]:
                    self.xyz[i, j] -= move[j]
                elif self.xyz[i, j] < limit_down[j]:
                    self.xyz[i, j] += move[j]
                else:
                    continue
        
    def write(self, locout):
        """
        Note: The cell origin may not be correct
        """
        
        with open(locout, 'w') as wt:
            wt.write('%d\n' % (self.TotNum))
            wt.write('Lattice=" %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f " '%(self.boxsize_x,0.0,0.0,0.0,self.boxsize_y,0.0,0.0,0.0,self.boxsize_z))
            wt.write('Origin=" %.5f %.5f %.5f "\n'%(-0.5*self.boxsize_x, -0.5*self.boxsize_y, -0.5*self.boxsize_z))
            for i in range(self.TotNum):
                wt.write('%s %.5f %.5f %.5f %s\n'%(self.element[self.elementSequence[i]], self.xyz[i,0],
                                                   self.xyz[i,1], self.xyz[i,2], self.elementSequence[i]))
        
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
        tmp=i.split('/')[-1]
        tmp=tmp.split('-')[3]
        file_number.append(int(tmp))
        
    file_number, file_name = zip(* sorted(zip(file_number, file_name) ) )
        
    print("The list of files:")
    
    for i in range(len(file_name)):
        print(file_number[i], file_name[i])
        
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
                    
                    if 'Origin' in linesplit[2]:
                        linesplit = linesplit[1].split()
                    else:
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

def mini_x_y(prog, ctype, file, eleCol, element, xDataCol, sizex, sizey, sizez):
    """
    Get full minimum positions of a cell by using the FCC analyzer.
    Arguments: the arguments of fccAnalyzer
    
    Return: 
        min_x, min_z: minimum positions of x and z planes;
        max_z: maximum position of z plane
        min_y: minimum position of y row on the minimum x plane;
        min_x_minZ: minimum position of x row on the minimum z plane;
        min_x_maxZ: minimum position of x row on the maximum z plane
    """
    screen = subprocess.run([prog, ctype, file, '%s'%eleCol, element,'%s'%xDataCol,\
                             '%s'%sizex, '%s'%sizey, '%s'%sizez],capture_output=True, text=True)
    
    screen = screen.stdout        
    lines = screen.split('\n')
    for i in range(len(lines)):
        
        if "Minimum lattice position" in lines[i]:
            min_x = float(lines[i].split()[3])
            min_y = float(lines[i].split()[4])
            min_z = float(lines[i].split()[5])
            max_z = float(lines[i].split()[6])
            
        if "minimum z plane:" in lines[i]:
            min_x_minZ = float(lines[i].split()[9])
            
        if "maximum z plane:" in lines[i]:
            min_x_maxZ = float(lines[i].split()[9])
            
        if "Reallocating" in lines[i]:
            continue
        elif "Allocating" in lines[i]:
            continue
        else:
            print(lines[i])
                
    return(min_x, min_y, min_z, max_z, min_x_minZ, min_x_maxZ)

def writeLmpIn(locout, name, x1, x2, y1, y2, z1, z2, time=1000):
    """
    Write an input file for LAMMPS.
    locout: output path for the input file
    name: the name of target cell without the .xyz format and the name for the dump file
    x1, x2, y1, y2, z1, z2: boundaries on three direactions [A]
    time: simulation total time step
    """
    
    with open(locout, "w") as wt:
        wt.write("variable \t name string '%s' \n" % name)
        wt.write("\n")
        wt.write("units \t metal \n")
        wt.write("atom_style	 \t atomic \n")
        wt.write("boundary \t p p p \n")
        wt.write("\n")
        wt.write("region	 \t mybox block %.5f %.5f %.5f %.5f %.5f %.5f units box \n" % \
                 (x1, x2, y1, y2, z1, z2))
        wt.write("create_box	 \t 1 mybox\n")
        wt.write("\n")
        wt.write("log \t log.${name} \n")
        wt.write("\n")
        wt.write("read_dump \t ${name}.xyz 0 x y z box no replace no trim no purge yes add yes format xyz \n")
        wt.write("\n")
        wt.write("pair_style	 \t eam/fs \n")
        if remote == "no":
            wt.write("pair_coeff	 \t * * /home/jinxinjx/Licument/Fusion_proA/W_potentials/W_AT-ZN_lammps.eam.fs W \n")
        elif remote == "yes":
            wt.write("pair_coeff	 \t * * /proj/jinxinjx/MD_Potentials/W_AT-ZN_lammps.eam.fs W \n")
        else:
            raise ValueError("the remote should be yes or no.")
        wt.write("\n")
        wt.write("thermo	 \t 100 \n")
        wt.write("thermo_style \t custom step temp press pe etotal dt time \n")
        wt.write("\n")
        wt.write("velocity \t all create 290 666 mom yes rot yes dist gaussian \n")
        wt.write("fix \t 1 all npt temp 290 290 0.1 iso 0.0 0.0 0.1 \n")
        wt.write("timestep \t 0.001 \n")
        wt.write("run \t %d\n" % time)
        wt.write("\n")
        wt.write("write_dump	 \t all custom ${name}.xyz element x y z type modify element W\n")
        
def readLmpCell(inpath, element):
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
                eleSequence.append(element)
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

def connect(path2FredericCellList, path2Base, prog, deltaZ, lmp_time, \
            nprocessor, lmp_prog, nprint):
    """
    The main function for connecting cells.
    Perform a relax after each connection
    
    path2FredericCellList: path to the folder of Frederic's cell
    path2Base: path to the folder of Base
    prog: the name of the fccAnalyzer program
    deltaZ: the margin distance between two cells
    lmp_time: the total time step used in LAMMPS simulations
    nprocessor: number of processors;
    lmp_prog: path to the LAMMPS program
    nprint: save the cell every nprint steps (useless)

    """

    time_c1 = time.time()
    
    #1. Get the list of Fredric's cells (sorted according to the order along the z axis)
    #   and get the average box sizes
    cellList = getFredricCellList(path2FredericCellList, character = '*.xyz.gz')
    prog = path2Base + prog
    
    x_list = []
    y_list = []
    
    for i in range(0, len(cellList)):
        dump1, tmp_x, tmp_y, tmp_z, dump2, dump3 = readFredric(cellList[i], 'W')
        x_list.append(tmp_x)
        y_list.append(tmp_y)
                
        #tmp_min_x, tmp_min_y = mini_x_y(prog, '-bcc', '-weakP', cellList[i], 5, 'W', 2, tmp_x, tmp_y, tmp_z)
        #mx_list.append(tmp_min_x)
        #my_list.append(tmp_min_y)
    
    x_list = np.array(x_list)
    y_list = np.array(y_list)    
    
    x_ave = np.average(x_list)
    y_ave = np.average(y_list)
            
            
    #2. Create two cells, at1 is the main cell, at2 will be connected to at1
    
    W_a_z = 2.73664028 #Lattice parameter of <111>-W along the z direction
    W_a_x = 7.74038759 #Lattice parameter of <111>-W along the x direction

    at1 = CELL()
    at2 = CELL()
    
    at1.element = {"W" : 1}
    at2.element = {"W" : 1}
    
    #3. Read the first Frederic's cell and extract all necessary information
    #   and expand the cell according to the average box x-y-sizes
    at1.TotNum, at1.boxsize_x, at1.boxsize_y, at1.boxsize_z, at1.xyz, at1.elementSequence = \
        readFredric(cellList[0], list(at1.element.keys())[0])
        
    at1.periodicBC(-at1.boxsize_x/2, at1.boxsize_x/2, -at1.boxsize_y/2, at1.boxsize_y/2, -at1.boxsize_z/2, at1.boxsize_z/2)
    
    #   (origin: corner)    
    at1.translate(at1.boxsize_x/2, at1.boxsize_y/2, 0)
    at1.shrink(x_ave/at1.boxsize_x, y_ave/at1.boxsize_y, 1)
    #   (origin: center) 
    at1.translate(-1*at1.boxsize_x/2, -1*at1.boxsize_y/2, 0) 
        
    #4. Get the minimum x and y latice positions of the first Frederic's cell by using the FCC analyzer
    #path2Base = '../../Base/'
    #prog = 'minixy'
    
    at1.write('tmp.xyz')
    at1.min_x, at1.min_y, at1.min_z, at1.max_z, at1.xz1, at1.xz2 = mini_x_y(prog, '-bcc', 'tmp.xyz', 5, 'W', 2, at1.boxsize_x, at1.boxsize_y, at1.boxsize_z)
    
    
    #Determine the type of the z plane (A, B or C),
    #according to the distance of the minimum x position on that plane with the globle minimum x position.
    #This applies to both the minimum z and maximum z planes
    row_distance1 = abs(abs(at1.xz1 -at1.min_x) % (W_a_x/2) )
    if row_distance1  > (W_a_x/4):
        head_z_mini_type = 'C'
        print("Type of minimum z plane: C")
    elif row_distance1 <= (W_a_x/4) and row_distance1 > (W_a_x/12):
        head_z_mini_type = 'B'
        print("Type of minimum z plane: B")
    else:
        head_z_mini_type = 'A'
        print("Type of minimum z plane: A")
    
    row_distance2 = abs(abs(at1.xz2 -at1.min_x) % (W_a_x/2) )    
    if row_distance2  > (W_a_x/4):
        head_z_max_type = 'C'
        print("Type of maximum z plane: C")
    elif row_distance2 <= (W_a_x/4) and row_distance2 > (W_a_x/12):
        head_z_max_type = 'B'
        print("Type of maximum z plane: B")
    else:
        head_z_max_type = 'A'
        print("Type of maximum z plane: A")   
    
    for i in range(1, len(cellList)):
        
        print("About to connect %s ...\n" % cellList[i])
        
        #5. Read the next Frederic's cell and extract all necessary information
        at2.TotNum, at2.boxsize_x, at2.boxsize_y, at2.boxsize_z, at2.xyz, at2.elementSequence = \
            readFredric(cellList[i], list(at1.element.keys())[0])
            
        #6. Get the minimum x and y latice positions of the next Frederic's cell
            #(First, Adjust the at2 boxsize on the x and y directions accroding to those of at1)
        at2.periodicBC(-at2.boxsize_x/2, at2.boxsize_x/2, -at2.boxsize_y/2, at2.boxsize_y/2, -at2.boxsize_z/2, at2.boxsize_z/2)  
        at2.translate(at2.boxsize_x/2, at2.boxsize_y/2, 0)
        at2.shrink((at1.boxsize_x/at2.boxsize_x), (at1.boxsize_y/at2.boxsize_y), 1)
        at2.translate(-1*at2.boxsize_x/2, -1*at2.boxsize_y/2, 0)    
            #(Then, use the fccAnalyzer)
        at2.write('tmp.xyz')
        at2.min_x, at2.min_y, at2.min_z, at2.max_z, at2.xz1, at2.xz2 = mini_x_y(prog, '-bcc', 'tmp.xyz', 5, 'W', 2, at2.boxsize_x, at2.boxsize_y, at2.boxsize_z)
        
        #7. Connect at1 cell (main cell) with the at2 cell
            #1). Move the origin of the at1 (main cell) to the corner (origin: corner)
        at1.translate(at1.boxsize_x/2, at1.boxsize_y/2, at1.boxsize_z/2)
        
            #2). Move the at2 cell on the xy plane according to the difference on the minimum lattice positions (origin: center)
        at2.translate((at1.min_x - at2.min_x), (at1.min_y - at2.min_y), 0)
        at2.periodicBC(-at2.boxsize_x/2, at2.boxsize_x/2, -at2.boxsize_y/2, at2.boxsize_y/2, -at2.boxsize_z/2, at2.boxsize_z/2)
        
        at2.write('tmp.xyz')
        at2.min_x, at2.min_y, at2.min_z, at2.max_z, at2.xz1, at2.xz2 = mini_x_y(prog, '-bcc', 'tmp.xyz', 5, 'W', 2, at2.boxsize_x, at2.boxsize_y, at2.boxsize_z)
        
            #Determine the type of plane on the at2 cell (tail),
            #The standard for this is the at1 cell (head)
        row_distance1 = abs(abs(at2.xz1 -at1.min_x) % (W_a_x/2) )
        if row_distance1  > (W_a_x/4):
            tail_z_mini_type = 'C'
            print("Type of minimum z plane: C")
        elif row_distance1 <= (W_a_x/4) and row_distance1 > (W_a_x/12):
            tail_z_mini_type = 'B'
            print("Type of minimum z plane: B")
        else:
            tail_z_mini_type = 'A'
            print("Type of minimum z plane: A")
        
        row_distance2 = abs(abs(at2.xz2 -at1.min_x) % (W_a_x/2) )    
        if row_distance2  > (W_a_x/4):
            tail_z_max_type = 'C'
            print("Type of maximum z plane: C")
        elif row_distance2 <= (W_a_x/4) and row_distance2 > (W_a_x/12):
            tail_z_max_type = 'B'
            print("Type of maximum z plane: B")
        else:
            tail_z_max_type = 'A'
            print("Type of maximum z plane: A")
            
            #If there is a stacking fault because of connection problem along the z direction,
            #we can move the at2 according to the following method.
            #Note in these cells from surface to bulk, the sequence of planes is like this:
            #CBA CBA CBA ...
            #This is important, since if the sequence is like ABC, then we need to use a different method
        if head_z_max_type == 'A' and tail_z_mini_type == 'B':
            at2.translate(0, 0, -1*(at2.min_z + at2.boxsize_z/2 + W_a_z/2)) 
        elif head_z_max_type == 'A' and tail_z_mini_type == 'A':
            at2.translate(0, 0, -1*(at2.min_z + at2.boxsize_z/2 + W_a_z/4)) 
        elif head_z_max_type == 'B' and tail_z_mini_type == 'C':
            at2.translate(0, 0, -1*(at2.min_z + at2.boxsize_z/2 + W_a_z/2)) 
        elif head_z_max_type == 'B' and tail_z_mini_type == 'B':
            at2.translate(0, 0, -1*(at2.min_z + at2.boxsize_z/2 + W_a_z/4)) 
        elif head_z_max_type == 'C' and tail_z_mini_type == 'A':
            at2.translate(0, 0, -1*(at2.min_z + at2.boxsize_z/2 + W_a_z/2))
        elif head_z_max_type == 'C' and tail_z_mini_type == 'C':
            at2.translate(0, 0, -1*(at2.min_z + at2.boxsize_z/2 + W_a_z/4)) 
            
        at2.periodicBC(-at2.boxsize_x/2, at2.boxsize_x/2, -at2.boxsize_y/2, at2.boxsize_y/2, -at2.boxsize_z/2, at2.boxsize_z/2)
        at2.write('tmp.xyz')
        
        
        at2.min_x, at2.min_y, at2.min_z, at2.max_z, at2.xz1, at2.xz2 = mini_x_y(prog, '-bcc', 'tmp.xyz', 5, 'W', 2, at2.boxsize_x, at2.boxsize_y, at2.boxsize_z)
        
        #Update the type of the plane 69792541
        #This will be the head cell for the next at.2 cell
        row_distance1 = abs(abs(at2.xz1 -at1.min_x) % (W_a_x/2) )
        if row_distance1  > (W_a_x/4):
            tail_z_mini_type = 'C'
            print("Type of minimum z plane: C")
        elif row_distance1 <= (W_a_x/4) and row_distance1 > (W_a_x/12):
            tail_z_mini_type = 'B'
            print("Type of minimum z plane: B")
        else:
            tail_z_mini_type = 'A'
            print("Type of minimum z plane: A")
            
        row_distance2 = abs(abs(at2.xz2 -at1.min_x) % (W_a_x/2) )   
        if row_distance2  > (W_a_x/4):
            head_z_max_type = 'C'
            print("Type of maximum z plane: C")
        elif row_distance2 <= (W_a_x/4) and row_distance2 > (W_a_x/12):
            head_z_max_type = 'B'
            print("Type of maximum z plane: B")
        else:
            head_z_max_type = 'A'
            print("Type of maximum z plane: A")
        
        
            #3). Move the at2 cell after the at1 cell (plus a small margin) (origin: x/y-center, z-ignore)
            #If the deltaZ is too small, there will be collision (indicated by the increase of temperature);
            #If the deltaZ is too big, MD code will spend efforts to reduce the volume
        #deltaZ = 0.5
        at2.translate(at2.boxsize_x/2, at2.boxsize_y/2, at2.boxsize_z/2 + at1.boxsize_z + deltaZ)
            #4). Connection (should be faster than numpy concantenate) (origin: center))
        at1.TotNum += at2.TotNum
        
        if i == (len(cellList) - 1):
            at1.boxsize_z += (at2.boxsize_z + deltaZ + deltaZ)
        else:
            at1.boxsize_z += (at2.boxsize_z + deltaZ)
        
        tmp_cell = np.empty([at1.TotNum, 3])
        tmp_cell[0 : at1.TotNum - at2.TotNum] = at1.xyz
        tmp_cell[at1.TotNum - at2.TotNum : ] = at2.xyz
        at1.xyz = tmp_cell
        tmp_cell = np.empty([1, 3])
        
        tmp_sequence = []
        #tmp_sequence[0 : at1.TotNum - at2.TotNum] = at1.elementSequence
        #tmp_sequence[at1.TotNum - at2.TotNum : ] = at2.elementSequence
        tmp_sequence.append(at1.elementSequence)
        tmp_sequence.append(at2.elementSequence)
        tmp_sequence = np.array(tmp_sequence)
        tmp_sequence = np.hstack(tmp_sequence)
        at1.elementSequence = tmp_sequence
        tmp_sequence = []
        
        at1.centerOrigin()
        
        at1.periodicBC(-at1.boxsize_x/2, at1.boxsize_x/2, -at1.boxsize_y/2, at1.boxsize_y/2, -at1.boxsize_z/2, at1.boxsize_z/2)
        
    
    at1.write('tmp.xyz')
    
    
    
    #8. Write the LAMMPS input file
    lmp_in = "lmp.in"
    lmp_out = "tmp"
    #lmp_time = 200
    writeLmpIn(lmp_in, lmp_out, -1*at1.boxsize_x/2, at1.boxsize_x/2, \
               -1*at1.boxsize_y/2, at1.boxsize_y/2, -1*at1.boxsize_z/2, at1.boxsize_z/2, lmp_time)
        
    #9. Relax the cell by using LAMMPS
    #nprocessor = 3
    #lmp_prog = '/home/jinxinjx/Licument/MD_Simu/lammps-29Oct20/src/lmp_mpi'
    subprocess.run(['mpirun', '--mca', 'btl_openib_allow_ib', \
                    '1', '-np', '%s'%nprocessor, lmp_prog, '-in', './%s'%lmp_in])#,capture_output=True, text=True)
        
    #10. Update the at1 cell by reading LAMMPS's dump file
        #1). Update the total number, coordinates and box sizes
    at1.TotNum, at1.boxsize_x, at1.boxsize_y, at1.boxsize_z, at1.xyz, at1.elementSequence = \
        readLmpCell('tmp.xyz', list(at1.element.keys())[0])
        #2). Move origin to the center (the origin of the LAMMPS output is very close the center)
    at1.periodicBC(-at1.boxsize_x/2, at1.boxsize_x/2, -at1.boxsize_y/2, at1.boxsize_y/2, -at1.boxsize_z/2, at1.boxsize_z/2)
        #3). Write the cell 
    at1.write('final.xyz')
    
    
    time_c2 = time.time()
    time_ch, time_cmin, time_cs = t2hms(time_c2 - time_c1)
    print("The connection of %s is finished (time: %d h %d min %.1f s).\n" % (cellList[i], time_ch, time_cmin, time_cs) )

#=============================================================================#

"""Initilization parameters"""

t1 = time.time()

path2FredericCellList = '../Project-3/AT-1-0.02dpa-KP/' 
path2Base = '../Base/'
prog = 'minixy_e'
deltaZ = 0.6
lmp_time = 50000
nprocessor = 50
if remote == 'no':
    lmp_prog = '/home/jinxinjx/Licument/MD_Simu/lammps-29Oct20/src/lmp_mpi'
elif remote == 'yes':
    lmp_prog = path2Base + 'lmp_mpi'
else:
    raise ValueError("the remote should be yes or no.")
nprint = 2

connect(path2FredericCellList, path2Base, prog, deltaZ, lmp_time, nprocessor, lmp_prog, nprint)

t2 = time.time()
th, tmin, ts = t2hms(t2 - t1)
print("Total time: %d h %d min %.1f s" % (th, tmin, ts))






