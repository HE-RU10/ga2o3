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
import peakutils
#import matplotlib.pyplot as plt

remote = "no"


class CELL:
    """
    """
    
    def __init__(self):
        
        self.element = {"H": 1}
        self.mass = {"H": 1}
        self.elementSequence = []
        self.massSeqeunce = []
        
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
                
    def peak(self, binWidth, peak_thres, peak_dist, figures=False):
        """
        Find the stacking peaks

        Parameters
        ----------
        binWidth : the width of a bin [A]
        peak_thres: threshold of peak search
        peak_dist: minimum distance between two bins

        Returns
        -------
        Minimum x plane
        Minimum y row
        Minimum z plane
        Maximum z plane
        Max x in min z 
        Max x in max z 

        """
        
        print("Performing peak analysis ...")
        
            
        binNum_x = int(self.boxsize_x / binWidth)
        binNum_y = int(self.boxsize_y / binWidth)
        binNum_z = int(self.boxsize_z / binWidth)
        
        binWidth_x = self.boxsize_x / binNum_x
        binWidth_y = self.boxsize_y / binNum_y
        binWidth_z = self.boxsize_z / binNum_z
        
        #1. Find peaks of x (planes)
        #bining the results
        counts_x, bins_x = np.histogram(self.xyz[:,0], binNum_x, (-self.boxsize_x/2, self.boxsize_x/2))
        peak_bins_x = peakutils.indexes(counts_x, peak_thres, peak_dist)
        peak_val_x = bins_x[peak_bins_x]
        #print(peak_val_x)
        
        min_x = peak_val_x[0]
        
        #2. Find peaks of y (rows in minimum x planes)
        #   2.1 inds gives the bin number of each coordinate
        inds = np.digitize(self.xyz[:,0], bins_x)
        
        #   Get all the coordinates on the minimum x planes
        val_min_x = []
        for i in range(len(inds)):
            if inds[i] == peak_bins_x[0]:
                val_min_x.append(self.xyz[i])
                
        val_min_x = np.array(val_min_x)
        # print(val_min_x)
        # print(len(val_min_x))
        
        #   2.2 Find the minimum y row
        counts_y, bins_y = np.histogram(val_min_x[:,1], binNum_y, (-self.boxsize_y/2, self.boxsize_y/2))
        # print(counts_y)
        # print(np.sum(counts_y))
        peak_bins_y = peakutils.indexes(counts_y, peak_thres, peak_dist)
        peak_val_y = bins_y[peak_bins_y]
        # print(peak_val_y)
        #print(bins_y)
        
        min_y = peak_val_y[0]
        
        #Plot the atoms on the minimum x plane
        if figures == True:
            fig1, axs = plt.subplots() 
            fig1.set_size_inches(4,4)
            axs.set_ylabel('Y',fontsize=14)
            axs.set_xlabel('Z',fontsize=14)
            axs.set_ylim(-self.boxsize_y/2, self.boxsize_y/2)
            axs.set_xlim(-self.boxsize_z/2, self.boxsize_z/2)
    
            axs.scatter(val_min_x[:,2], val_min_x[:,1], s=10)
            for i in range(len(peak_val_y)):
                axs.plot([np.min(val_min_x[:,2]), np.max(val_min_x[:,2])], [peak_val_y[i], peak_val_y[i]], color='C1', linewidth=1)
    
            axs.grid(visible=True, which='major',linestyle='--')
            plt.show()
            fig1.tight_layout()
            fig1.savefig('mini_x_plane.png', dpi=300)#,facecolor='none')
            
            
        #3. Find peaks of z (planes)
        #bining the results
        #Add one more bin by adding "+1" and "+binWidth_z", because peakutils cannot recongnize the peak in the last bin
        counts_z, bins_z = np.histogram(self.xyz[:,2], binNum_z+1, (-self.boxsize_z/2, self.boxsize_z/2+binWidth_z))
        #print(counts_z)
        peak_bins_z = peakutils.indexes(counts_z, peak_thres, peak_dist)
        peak_val_z = bins_z[peak_bins_z]
        peak_counts_z = counts_z[peak_bins_z]
        #print(peak_bins_z)
        #print(peak_val_z)
        
        
        min_z = peak_val_z[0]
        max_z = peak_val_z[-1]
        
        #Plot the peaks of z planes
        if figures == True:
            fig1, axs = plt.subplots() #, gridspec_kw = {'hspace' : 0, 'height_ratios': [2, 2, 0.7, 1.5]})
            fig1.set_size_inches(8,4)
            axs.set_ylabel('Counts',fontsize=14)
            axs.set_xlabel('Z',fontsize=14)

            axs.hist(bins_z[:-1], bins_z, weights=counts_z, histtype='bar', label='Inversion', color='C2', edgecolor='k') #In [-10,10]
            axs.scatter(peak_val_z, peak_counts_z, s=4)
            plt.show()
            fig1.tight_layout()
            fig1.savefig('peaks_z_planes', dpi=300)#,facecolor='none')

        
        #4. Find peaks of x (rows on minimum and maximum z planes)
        #   4.1 inds gives the bin number of each coordinate
        inds = np.digitize(self.xyz[:,2], bins_z)
        
        #   Get all the coordinates on the minimum z planes
        val_min_z = []
        val_max_z = []
        for i in range(len(inds)):
            if inds[i] == peak_bins_z[0]:
                val_min_z.append(self.xyz[i])
            if inds[i] == peak_bins_z[-1]:
                val_max_z.append(self.xyz[i])
                
        val_min_z = np.array(val_min_z)
        val_max_z = np.array(val_max_z)
        
        #print(val_min_z)
        # print(val_max_z)
        
        #   4.2 Find the minimum x row on the minimum and maximum z planes
        counts_x_minz, bins_x_minz = np.histogram(val_min_z[:,0], binNum_x+1, (-self.boxsize_x/2, self.boxsize_x/2+binWidth_x))
        counts_x_maxz, bins_x_maxz = np.histogram(val_max_z[:,0], binNum_x+1, (-self.boxsize_x/2, self.boxsize_x/2)+binWidth_x)
        #print(counts_x_maxz)

        #Set the threshold to 0.99, so only take the maximum peak
        peak_bins_x_minz = peakutils.indexes(counts_x_minz, 0.99, peak_dist)
        peak_bins_x_maxz = peakutils.indexes(counts_x_maxz, 0.99, peak_dist)
        peak_val_x_minz = bins_x_minz[peak_bins_x_minz]
        peak_val_x_maxz = bins_x_maxz[peak_bins_x_maxz]
        
        # print(peak_val_x_minz)
        # print(peak_val_x_maxz)
        # print(peak_bins_x_maxz)
        # print(counts_x_maxz[peak_bins_x_maxz])
        
        maxPeak_val_x_minz = peak_val_x_minz[0]
        maxPeak_val_x_maxz = peak_val_x_maxz[0]
        
        #Plot the atoms on the minimum z plane
        if figures == True:
            fig1, axs = plt.subplots() 
            fig1.set_size_inches(4,4)
            axs.set_ylabel('Y',fontsize=14)
            axs.set_xlabel('X',fontsize=14)
            axs.set_ylim(-self.boxsize_y/2, self.boxsize_y/2)
            axs.set_xlim(-self.boxsize_x/2, self.boxsize_x/2+binWidth_x)
    
            axs.scatter(val_min_z[:,0], val_min_z[:,1], s=10)
            for i in range(len(peak_val_x_minz)):
                axs.plot([peak_val_x_minz[i], peak_val_x_minz[i]], [-self.boxsize_y/2, self.boxsize_y/2], color='C1', linewidth=1)
    
            axs.grid(visible=True, which='major',linestyle='--')
            plt.show()
            fig1.tight_layout()
            fig1.savefig('mini_z_plane.png', dpi=300)#,facecolor='none')
            
            
        #Plot the atoms on the maximum z plane
        if figures == True:
            fig1, axs = plt.subplots() 
            fig1.set_size_inches(4,4)
            axs.set_ylabel('Y',fontsize=14)
            axs.set_xlabel('X',fontsize=14)
            axs.set_ylim(-self.boxsize_y/2, self.boxsize_y/2)
            axs.set_xlim(-self.boxsize_x/2, self.boxsize_x/2+binWidth_x)
    
            axs.scatter(val_max_z[:,0], val_max_z[:,1], s=10)
            for i in range(len(peak_val_x_maxz)):
                axs.plot([peak_val_x_maxz[i], peak_val_x_maxz[i]], [-self.boxsize_y/2, self.boxsize_y/2], color='C1', linewidth=1)
    
            axs.grid(visible=True, which='major',linestyle='--')
            plt.show()
            fig1.tight_layout()
            fig1.savefig('max_z_plane.png', dpi=300)#,facecolor='none')


        print('Minimum x plane: %7.3f' % min_x)
        print('Minimum y row  : %7.3f' % min_y)
        print('Minimum z plane: %7.3f' % min_z)
        print('Maximum z plane: %7.3f' % max_z)
        print('Max x in min z : %7.3f' % maxPeak_val_x_minz)
        print('Max x in max z : %7.3f' % maxPeak_val_x_maxz)
       
        return(min_x, min_y, min_z, max_z, maxPeak_val_x_minz, maxPeak_val_x_maxz)
    
    def hist(self, counts, bins, locout='tmp.png'):
        """
        Plot histograms

        """
        
        fig1, axs = plt.subplots() #, gridspec_kw = {'hspace' : 0, 'height_ratios': [2, 2, 0.7, 1.5]})
        #fig1=plt.figure()
        fig1.set_size_inches(8,4)
        #axs1=fig1.add_subplot(111)
        axs.set_ylabel('Counts',fontsize=14)
        #axs1.set_ylim(0, 9)


        #axs[0].hist(res_inv_bins[:-1], res_inv_bins, weights=res_inv_counts, histtype='bar', label='Inversion', color='C2', edgecolor='k') #In [-inf,inf]
        axs.hist(bins[:-1], bins, weights=counts, histtype='bar', label='Inversion', color='C2', edgecolor='k') #In [-10,10]
        #axs[0].plot(x, Fx/5, label='F(x)/5', color='C1') #In [-inf,inf]
        plt.show()
        fig1.tight_layout()
        #fig1.savefig('tmp.png',dpi=600)#,facecolor='none')
        fig1.savefig('%s'%locout, dpi=400)#,facecolor='none')
        
        
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
                
    def write_lmp(self, locout):
        """
        Write output that can be read by LAMMPS read_data command

        """
        
        print("\nWriting output for LAMMPS: %s\n" % locout)
        
        with open(locout, 'w') as op:
            op.write("# LAMMPS data file written by TLG\n")
            op.write("%d atoms\n" % (self.TotNum))
            op.write("%d atom types\n" % (len(self.element)))
            op.write("%.5f %.5f xlo xhi\n" % (-0.5*self.boxsize_x, 0.5*self.boxsize_x))
            op.write("%.5f %.5f ylo yhi\n" % (-0.5*self.boxsize_y, 0.5*self.boxsize_y))
            op.write("%.5f %.5f zlo zhi\n" % (-0.5*self.boxsize_z, 0.5*self.boxsize_z))
            op.write("\n")
            op.write("Masses\n")
            op.write("\n")
            for i in range(len(self.mass)):
                op.write("%d %.5f   #%s\n" % (i+1, self.mass[self.massSeqeunce[i]], self.massSeqeunce[i]))
            op.write("\n")
            op.write("Atoms\n")
            op.write("\n")
            for i in range(len(self.xyz)):
                op.write("%d %d %.5f %.5f %.5f\n" % (i+1, self.element[self.elementSequence[i]], self.xyz[i,0], self.xyz[i,1], self.xyz[i,2]))
                
    def readSelf(self, locin):
        """
        Read the cell written by the function write()
        """
        
        coord = []
        elementSequence = []
        
        with open(locin, 'r') as op:
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
                    
        self.TotNum = TotNum
        self.boxsize_x = boxsize_x
        self.boxsize_y = boxsize_y
        self.boxsize_z = boxsize_z
        self.elementSequence = np.array(elementSequence)
        self.xyz = np.array(coord)
        
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

def getHeRuCellList(path, character = '*/data.*'):
    """
    Get the list of He Ru's cells and sorted according to the number along the z axis.
    
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
        tmp=tmp.split('-')[-1]
        file_number.append(int(tmp))
        
    file_number, file_name = zip(* sorted(zip(file_number, file_name) ) )
        
    print("The list of files:")
    
    for i in range(len(file_name)):
        print(file_number[i], file_name[i])
        
    return(file_name)

def readASE(inpath, elementID=None):
    """
    Read the cells generated by ASE
    inpath: path to the cell;
    elementID: if not none,
               only take atoms in the list of element IDs
    """
    
    print("Start to read ASE cell: %s\n" % inpath)
    
    coord = []
    eleSequence = []    
    
    if '.gz' in inpath:
        
        with gzip.open(inpath, 'rt') as op:
            
            lineCount = 0
            line_flag = 1e10
            
            for line in op:
                lineCount += 1
                
                #line_flag is the line number where 'Atoms' is, and is followed by coordinates
                if 'Atoms' in line:
                    line_flag = lineCount
                
                if lineCount == 3:
                    totNum = int(line.split()[0])
                    
                if lineCount > 3 and 'Atoms' not in line:
                    if 'xlo' in line:
                        linesplit = line.split()
                        boxsize_x = float(linesplit[1]) - float(linesplit[0])
                    if 'ylo' in line:
                         linesplit = line.split()
                         boxsize_y = float(linesplit[1]) - float(linesplit[0])
                    if 'zlo' in line:
                         linesplit = line.split()
                         boxsize_z = float(linesplit[1]) - float(linesplit[0])  
                         
                if lineCount >= (line_flag+2) and lineCount < (line_flag+2+totNum):
                    linesplit = line.split()
                    if elementID == None:
                        eleSequence.append(int(linesplit[1]))
                        coord.append([])
                        coord[-1].append(float(linesplit[2]))
                        coord[-1].append(float(linesplit[3]))
                        coord[-1].append(float(linesplit[4]))
                    else:
                        if int(linesplit[1]) in elementID:
                            eleSequence.append(int(linesplit[1]))
                            coord.append([])
                            coord[-1].append(float(linesplit[2]))
                            coord[-1].append(float(linesplit[3]))
                            coord[-1].append(float(linesplit[4]))
                    
                if 'Velocities' in line:
                    break
                
    else:
        
        with open(inpath, 'r') as op:
            
            lineCount = 0
            line_flag = 1e10
            
            for line in op:
                lineCount += 1
                
                #line_flag is the line number where 'Atoms' is, and is followed by coordinates
                if 'Atoms' in line:
                    line_flag = lineCount
                
                if lineCount == 3:
                    totNum = int(line.split()[0])
                    
                if lineCount > 3 and 'Atoms' not in line:
                    if 'xlo' in line:
                        linesplit = line.split()
                        boxsize_x = float(linesplit[1]) - float(linesplit[0])
                    if 'ylo' in line:
                         linesplit = line.split()
                         boxsize_y = float(linesplit[1]) - float(linesplit[0])
                    if 'zlo' in line:
                         linesplit = line.split()
                         boxsize_z = float(linesplit[1]) - float(linesplit[0])  
                         
                if lineCount >= (line_flag+2) and lineCount < (line_flag+2+totNum):
                    linesplit = line.split()
                    if elementID == None:
                        eleSequence.append(int(linesplit[1]))
                        coord.append([])
                        coord[-1].append(float(linesplit[2]))
                        coord[-1].append(float(linesplit[3]))
                        coord[-1].append(float(linesplit[4]))
                    else:
                        if int(linesplit[1]) in elementID:
                            eleSequence.append(int(linesplit[1]))
                            coord.append([])
                            coord[-1].append(float(linesplit[2]))
                            coord[-1].append(float(linesplit[3]))
                            coord[-1].append(float(linesplit[4]))
                    
                if 'Velocities' in line:
                    break
     
    coord = np.array(coord)
    eleSequence = np.array(eleSequence)
    
    if elementID == None:
        if len(coord) != totNum: 
            raise ValueError("The total number of atoms does not match the number of coordinates")
    else:
        totNum = len(coord)
    
    return(totNum, boxsize_x, boxsize_y, boxsize_z, coord, eleSequence)

def mini_x_y(prog, ctype, file, eleCol, element, xDataCol, sizex, sizey, sizez, details=False):
    """
    Get full minimum positions of a cell by using the FCC analyzer.
    Arguments: the arguments of fccAnalyzer
    
    If details = Flase, do not print FCCanalyzer details
    else, print FCCanalyzer details
    
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
            if details == True:
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
            wt.write("pair_coeff	 \t * * eam/fs /home/heruhe/Desktop/Ga2o3/potential/Ga-O_20220921.eam.fs Ga Ga O O O \n")
            wt.write("pair_coeff	 \t * * tabgap /home/heruhe/Desktop/Ga2o3/potential/Ga-O_20220921.tabgap Ga Ga O O O no yes \n")
        elif remote == "yes":
            wt.write("pair_coeff	 \t * * eam/fs /proj/jinxinjx/MD_Potentials/Ga-O_20220921.eam.fs Ga Ga O O O \n")
            wt.write("pair_coeff	 \t * * tabgap /proj/jinxinjx/MD_Potentials/Ga-O_20220921.tabgap Ga Ga O O O no yes \n")
        else:
            raise ValueError("the remote should be yes or no.")
        wt.write("\n")
        wt.write("thermo	 \t 100 \n")
        wt.write("thermo_style \t custom step temp press pe etotal dt time \n")
        wt.write("\n")
        wt.write("min_style \t cg \n")
        wt.write("minimize \t 1e-5 1e-6 10000 10000 \n")
        wt.write("\n")
        wt.write("velocity \t all create 290 666 mom yes rot yes dist gaussian \n")
        wt.write("fix \t 1 all npt temp 290 290 0.1 iso 0.0 0.0 0.1 \n")
        wt.write("timestep \t 0.001 \n")
        wt.write("run \t %d\n" % time)
        wt.write("\n")
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

def readXYZ(inpath, elementSelect=None):
    """
    Read the classical .xyz cell
    
    elementSelect: if not none,
               only take atoms in the list of element IDs
               
               
    Return: cells with origins at zero

    """
    #1. Read the damaged cell
    coord = []
    elementSequence = []
    elementIDsequence = []
    
    with open(inpath, 'r') as op:
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
                origin_x = float(linesplit[12])
                origin_y = float(linesplit[13])
                origin_z = float(linesplit[14])
            else:
                if elementSelect == None:
                    elementSequence.append(linesplit[4])
                    elementIDsequence.append(int(linesplit[0]))
                    coord.append([])
                    coord[-1].append(float(linesplit[1]))
                    coord[-1].append(float(linesplit[2]))
                    coord[-1].append(float(linesplit[3]))
                else:
                    if linesplit[4] in elementSelect:
                        elementSequence.append(linesplit[4])
                        elementIDsequence.append(int(linesplit[0]))
                        coord.append([])
                        coord[-1].append(float(linesplit[1]))
                        coord[-1].append(float(linesplit[2]))
                        coord[-1].append(float(linesplit[3]))
                
    coord = np.array(coord)
    elementSequence = np.array(elementSequence)
    elementIDsequence = np.array(elementIDsequence)
    
    TotNum = len(coord)
    
    #2. Move the cell origin to zero
    coord -= np.array([origin_x, origin_y, origin_z])
    
    
    return(TotNum, boxsize_x, boxsize_y, boxsize_z, coord, elementSequence, elementIDsequence)
    
        
def t2hms(t):
    """
    Convert time t in seconds to hour-min-s
    """
    time_h = int ( t / 3600 )
    time_min = int (t / 60) - time_h * 60
    time_s = (t) % 60
    
    return(time_h, time_min, time_s)

#=============================================================================#

def connect(path2CellList, path2Base, prog, deltaZ, lmp_time, \
            nprocessor, lmp_prog, nprint, a_z, n_x):
    """
    The main function for connecting cells.
    Perform a relax after each connection
    
    path2CellList: path to the folder of MD cells
    path2Base: path to the folder of Base
    prog: the name of the fccAnalyzer program
    deltaZ: the margin distance between two cells
    lmp_time: the total time step used in LAMMPS simulations
    nprocessor: number of processors;
    lmp_prog: path to the LAMMPS program
    nprint: save the cell every nprint steps (useless)
    a_z: Lattice parameter along the z direction
    n_x: number of unit cell along the x direction

    """

    
    time_c1 = time.time()
    
    #1. Get the list of He Ru's cells (sorted according to the order along the z axis)
    #   and get the averag 
    cellList = getHeRuCellList(path2CellList, character = '*/data.*')
    prog = path2Base + prog
    
    
    x_list = []
    y_list = []
    
    for i in range(0, len(cellList)):
        dump1, tmp_x, tmp_y, tmp_z, dump2, dump3, dump4 = readXYZ(cellList[i])
        x_list.append(tmp_x)
        y_list.append(tmp_y)
                    
    x_list = np.array(x_list)
    y_list = np.array(y_list)    
    
    x_ave = np.average(x_list)
    y_ave = np.average(y_list)
                    
    a_x = x_ave / n_x
    
    #2. Create two cells, at1 is the main cell, at2 will be connected to at1

    #at1 and at2 are used by FCCanalyzer
    at1 = CELL()
    at2 = CELL()
    #at_con1, at_con2 are cells to be connected
    at_con1 = CELL()
    at_con2 = CELL()
    
    at1.element = {"Ga1" : 1, "Ga2" : 2, "O1" : 3, "O2" : 4, "O3" : 5}
    at2.element = {"Ga1" : 1, "Ga2" : 2, "O1" : 3, "O2" : 4, "O3" : 5}
    at_con1.element = {"Ga1" : 1, "Ga2" : 2, "O1" : 3, "O2" : 4, "O3" : 5}
    at_con2.element = {"Ga1" : 1, "Ga2" : 2, "O1" : 3, "O2" : 4, "O3" : 5}
    
    at_con1.mass = {"Ga1" : 69.723, "Ga2" : 69.723, "O1" : 15.9994, "O2" : 15.9994, "O3" : 15.9994}
    at_con1.massSeqeunce = ["Ga1", "Ga2", "O1", "O2", "O3"]
    
    
    #3. Read the first He Ru's cell and extract all necessary information
    #   !!!          only take the 3 Oxygen atoms
    #   and expand the cell according to the average box x-y-sizes
    at1.TotNum, at1.boxsize_x, at1.boxsize_y, at1.boxsize_z, at1.xyz, at1.elementSequence, at1.elementID = \
        readXYZ(cellList[0], ['O1', 'O2', 'O3'])  #For FCCanalyzer
        
    at_con1.TotNum, at_con1.boxsize_x, at_con1.boxsize_y, at_con1.boxsize_z, at_con1.xyz, at_con1.elementSequence, at_con1.elementID = \
        readXYZ(cellList[0])  #For connection
        
    at1.elementSequence = np.full((at1.TotNum), 'O')
    at1.element = {"O" : 1}
    
    #at1.periodicBC(-at1.boxsize_x/2, at1.boxsize_x/2, -at1.boxsize_y/2, at1.boxsize_y/2, -at1.boxsize_z/2, at1.boxsize_z/2)
    
    #   (origin: corner)    
    #at1.translate(at1.boxsize_x/2, at1.boxsize_y/2, 0)
    at1.shrink(x_ave/at1.boxsize_x, y_ave/at1.boxsize_y, 1)
    at_con1.shrink(x_ave/at_con1.boxsize_x, y_ave/at_con1.boxsize_y, 1)
    #   (origin: center) 
    at1.translate(-1*at1.boxsize_x/2, -1*at1.boxsize_y/2, -1*at1.boxsize_z/2) 
    at_con1.translate(-1*at_con1.boxsize_x/2, -1*at_con1.boxsize_y/2, -1*at_con1.boxsize_z/2) 
        
    #4. Get the minimum x and y latice positions of the first cell by using the FCC analyzer
    #path2Base = '../../Base/'
    #prog = 'minixy'
    
    at1.write('tmp.xyz')
    
    print('\nGet the minimum lattice positions of 1st cell ...')
    #at1.min_x, at1.min_y, at1.min_z, at1.max_z, at1.xz1, at1.xz2 = mini_x_y(prog, '-bcc', 'tmp.xyz', 5, 'O', 2, at1.boxsize_x, at1.boxsize_y, at1.boxsize_z, details=True)
    
    peak_thres = 0.5
    peak_dist = 0.5
    peak_bin = 0.2
    at1.min_x, at1.min_y, at1.min_z, at1.max_z, at1.xz1, at1.xz2 = at1.peak(peak_bin, peak_thres, peak_dist, figures=False)

    
    #Determine the type of the z plane (A or B),
    #according to the distance of the minimum x position on that plane with the globle minimum x position.
    #This applies to both the minimum z and maximum z planes
    print("\nStacking sequence of first cell:")
    row_distance1 = abs(abs(at1.xz1 -at1.min_x) % a_x )
    if row_distance1 <= (a_x*3/4) and row_distance1 > (a_x/4):
        head_z_mini_type = 'B'
        print("Type of minimum z plane: B")
    else:
        head_z_mini_type = 'A'
        print("Type of minimum z plane: A")
    
    row_distance2 = abs(abs(at1.xz2 -at1.min_x) % a_x )    
    if row_distance2 <= (a_x*3/4) and row_distance2 > (a_x/4):
        head_z_max_type = 'B'
        print("Type of maximum z plane: B\n")
    else:
        head_z_max_type = 'A'
        print("Type of maximum z plane: A\n")       
    
    
    for i in range(1, len(cellList)):
        
        print("\nAbout to connect %s ..." % cellList[i])
        
        #5. Read the next cell and extract all necessary information
        at2.TotNum, at2.boxsize_x, at2.boxsize_y, at2.boxsize_z, at2.xyz, at2.elementSequence, at2.elementID = \
            readXYZ(cellList[i], ['O1', 'O2', 'O3'])   #For FCCanalyzer
            
        at_con2.TotNum, at_con2.boxsize_x, at_con2.boxsize_y, at_con2.boxsize_z, at_con2.xyz, at_con2.elementSequence, at_con2.elementID = \
            readXYZ(cellList[i])   #For connection
        
        at2.elementSequence = np.full((at2.TotNum), 'O')
        at2.element = {"O" : 1}
        
        #6. Get the minimum x and y latice positions of the next cell
            #(First, Adjust the at2 boxsize on the x and y directions accroding to those of at1)
        #at2.periodicBC(-at2.boxsize_x/2, at2.boxsize_x/2, -at2.boxsize_y/2, at2.boxsize_y/2, -at2.boxsize_z/2, at2.boxsize_z/2)  
        #at2.translate(at2.boxsize_x/2, at2.boxsize_y/2, 0)
        at2.shrink((at1.boxsize_x/at2.boxsize_x), (at1.boxsize_y/at2.boxsize_y), 1)
        at_con2.shrink((at_con1.boxsize_x/at_con2.boxsize_x), (at_con1.boxsize_y/at_con2.boxsize_y), 1)
        
        at2.translate(-1*at2.boxsize_x/2, -1*at2.boxsize_y/2, -1*at2.boxsize_z/2)    
        at_con2.translate(-1*at_con2.boxsize_x/2, -1*at_con2.boxsize_y/2, -1*at_con2.boxsize_z/2)    
            #(Then, use the fccAnalyzer)
        at2.write('tmp.xyz')
        print('\nAdjusting lattice position ...')
        #at2.min_x, at2.min_y, at2.min_z, at2.max_z, at2.xz1, at2.xz2 = mini_x_y(prog, '-bcc', 'tmp.xyz', 5, 'O', 2, at2.boxsize_x, at2.boxsize_y, at2.boxsize_z, details=True)
        
        at2.min_x, at2.min_y, at2.min_z, at2.max_z, at2.xz1, at2.xz2 = at2.peak(peak_bin, peak_thres, peak_dist, figures=False)
        
        
        #7. Connect at1 cell (main cell) with the at2 cell
            #1). Move the origin of the at1 (main cell) to the corner (origin: corner)
        at_con1.translate(at_con1.boxsize_x/2, at_con1.boxsize_y/2, at_con1.boxsize_z/2)
        
            #2). Move the at2 cell on the xy plane according to the difference on the minimum lattice positions (origin: center)
        
        print('\nTranslate: %.5f %.5f\n' % ((at1.min_x - at2.min_x), (at1.min_y - at2.min_y)))
        at2.translate((at1.min_x - at2.min_x), (at1.min_y - at2.min_y), 0)
        at2.periodicBC(-at2.boxsize_x/2, at2.boxsize_x/2, -at2.boxsize_y/2, at2.boxsize_y/2, -at2.boxsize_z/2, at2.boxsize_z/2)
        
        at_con2.translate((at1.min_x - at2.min_x), (at1.min_y - at2.min_y), 0)
        at_con2.periodicBC(-at_con2.boxsize_x/2, at_con2.boxsize_x/2, -at_con2.boxsize_y/2, at_con2.boxsize_y/2, -at_con2.boxsize_z/2, at_con2.boxsize_z/2)
        
        at2.write('tmp.xyz')
        print('Checking stacking sequence (1) ...')
        #at2.min_x, at2.min_y, at2.min_z, at2.max_z, at2.xz1, at2.xz2 = mini_x_y(prog, '-bcc', 'tmp.xyz', 5, 'O', 2, at2.boxsize_x, at2.boxsize_y, at2.boxsize_z, details=True)
        
        at2.min_x, at2.min_y, at2.min_z, at2.max_z, at2.xz1, at2.xz2 = at2.peak(peak_bin, peak_thres, peak_dist, figures=False)
        
            #Determine the type of plane on the at2 cell (tail),
            #The standard for this is the at1 cell (head)
        row_distance1 = abs(abs(at2.xz1 -at1.min_x) % a_x )
        print()
        if row_distance1 <= (a_x*3/4) and row_distance1 > (a_x/4):
            tail_z_mini_type = 'B'
            print("Type of minimum z plane: B")
        else:
            tail_z_mini_type = 'A'
            print("Type of minimum z plane: A")
        
        row_distance2 = abs(abs(at2.xz2 -at1.min_x) % a_x )    
        if row_distance2 <= (a_x*3/4) and row_distance2 > (a_x/4):
            tail_z_max_type = 'B'
            print("Type of maximum z plane: B")
        else:
            tail_z_max_type = 'A'
            print("Type of maximum z plane: A")   
            
            #If there is a stacking fault because of connection problem along the z direction,
            #we can move the at2 according to the following method.
            #Note in these cells from surface to bulk, the sequence of planes is like this:
            #AB AB AB ...
            
        stackFlag = 1
            
        if head_z_max_type == 'A' and tail_z_mini_type == 'A':
            stackFlag = 0
            print('\n    Found a stacking fault AA')
            at2.translate(0, 0, -1*(at2.min_z + at2.boxsize_z/2 + a_z/4)) 
            at_con2.translate(0, 0, -1*(at2.min_z + at2.boxsize_z/2 + a_z/4)) 
            print('    Trying to solve it...\n')
        elif head_z_max_type == 'B' and tail_z_mini_type == 'B':
            stackFlag = 0
            print('\n    Found a stacking fault BB')
            at2.translate(0, 0, -1*(at2.min_z + at2.boxsize_z/2 + a_z/4)) 
            at_con2.translate(0, 0, -1*(at2.min_z + at2.boxsize_z/2 + a_z/4)) 
            print('    Trying to solve it...\n')
        
        at2.periodicBC(-at2.boxsize_x/2, at2.boxsize_x/2, -at2.boxsize_y/2, at2.boxsize_y/2, -at2.boxsize_z/2, at2.boxsize_z/2)
        at_con2.periodicBC(-at_con2.boxsize_x/2, at_con2.boxsize_x/2, -at_con2.boxsize_y/2, at_con2.boxsize_y/2, -at_con2.boxsize_z/2, at_con2.boxsize_z/2)
        
        at2.write('tmp.xyz')
        
        if stackFlag == 0:
            print('\nChecking stacking sequence (2) ...')
            #at2.min_x, at2.min_y, at2.min_z, at2.max_z, at2.xz1, at2.xz2 = mini_x_y(prog, '-bcc', 'tmp.xyz', 5, 'O', 2, at2.boxsize_x, at2.boxsize_y, at2.boxsize_z, details=True)
            at2.min_x, at2.min_y, at2.min_z, at2.max_z, at2.xz1, at2.xz2 = at2.peak(peak_bin, peak_thres, peak_dist, figures=False)
    
            #Update the type of the plane 
            #This will be the head cell for the next at.2 cell
            row_distance1 = abs(abs(at2.xz1 - at1.min_x) % a_x )
            print()
            if row_distance1 <= (a_x*3/4) and row_distance1 > (a_x/4):
                tail_z_mini_type = 'B'
                print("Type of minimum z plane (2nd check): B")
            else:
                tail_z_mini_type = 'A'
                print("Type of minimum z plane (2nd check): A")
                
            row_distance2 = abs(abs(at2.xz2 -at1.min_x) % a_x )   
            if row_distance2 <= (a_x*3/4) and row_distance2 > (a_x/4):
                head_z_max_type = 'B'
                print("Type of maximum z plane (2nd check): B")
            else:
                head_z_max_type = 'A'
                print("Type of maximum z plane (2nd check): A")
                
            stackFlag = 1
         

            #3). Move the at2 cell after the at1 cell (plus a small margin) (origin: x/y-center, z-ignore)
            #If the deltaZ is too small, there will be collision (indicated by the increase of temperature);
            #If the deltaZ is too big, MD code will spend efforts to reduce the volume
        at_con2.translate(at_con2.boxsize_x/2, at_con2.boxsize_y/2, at_con2.boxsize_z/2 + at_con1.boxsize_z + deltaZ)
            #4). Connection (should be faster than numpy concantenate) (origin: center))
        at_con1.TotNum += at_con2.TotNum
        
        if i == (len(cellList) - 1):
            at_con1.boxsize_z += (at_con2.boxsize_z + deltaZ + deltaZ)
        else:
            at_con1.boxsize_z += (at_con2.boxsize_z + deltaZ)
        
        tmp_cell = np.empty([at_con1.TotNum, 3])
        tmp_cell[0 : at_con1.TotNum - at_con2.TotNum] = at_con1.xyz
        tmp_cell[at_con1.TotNum - at_con2.TotNum : ] = at_con2.xyz
        at_con1.xyz = tmp_cell
        tmp_cell = np.empty([1, 3])
        
        tmp_sequence = []
        #tmp_sequence[0 : at1.TotNum - at2.TotNum] = at1.elementSequence
        #tmp_sequence[at1.TotNum - at2.TotNum : ] = at2.elementSequence
        for j in range(len(at_con1.elementSequence)):
            tmp_sequence.append(at_con1.elementSequence[j])
        for j in range(len(at_con2.elementSequence)):
            tmp_sequence.append(at_con2.elementSequence[j])
        tmp_sequence = np.array(tmp_sequence)
        at_con1.elementSequence = tmp_sequence
        tmp_sequence = []
        
        at_con1.centerOrigin()
        
        at_con1.periodicBC(-at_con1.boxsize_x/2, at_con1.boxsize_x/2, -at_con1.boxsize_y/2, at_con1.boxsize_y/2, -at_con1.boxsize_z/2, at_con1.boxsize_z/2)
        
         
    at_con1.write_lmp('tmp.xyz')
    
    
    
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
    at_con1.TotNum, at_con1.boxsize_x, at_con1.boxsize_y, at_con1.boxsize_z, at_con1.xyz, at_con1.elementSequence = \
        readLmpCell('tmp.xyz')
        #2). Move origin to the center (the origin of the LAMMPS output is very close to the center)
    at_con1.periodicBC(-at_con1.boxsize_x/2, at_con1.boxsize_x/2, -at_con1.boxsize_y/2, at_con1.boxsize_y/2, -at_con1.boxsize_z/2, at_con1.boxsize_z/2)
        #3). Write the cell 
    at_con1.write('final.xyz')
    
    
    time_c2 = time.time()
    time_ch, time_cmin, time_cs = t2hms(time_c2 - time_c1)
    print("The connection of %s is finished (time: %d h %d min %.1f s).\n" % (cellList[i], time_ch, time_cmin, time_cs) )
    

#=============================================================================#

"""Initilization parameters"""

if __name__ == "__main__":

    print("Runing merRelax.py")

    #!!!Please pay attention to the path under "remote" (some paths are in the functions)
    
    t1 = time.time()
    
    a_z = 3.094884375 #Lattice parameter of O sub-lattice along the z direction
    n_z = 32 #Number of unit cell along the z direction
    a_x = 4.047675    #Lattice parameter of O sub-lattice along the x direction
    n_x = 24 #Number of unit cell along the x direction
    
    if remote == 'no':
        path2HeRuCellList = '/home/heruhe/Desktop/Ga2o3/code/ga2o3/stacking_box/prepare_target/Rotate' 
        path2Base = '/home/heruhe/Desktop/Ga2o3/code/ga2o3/stacking_box/simu/Base/'
    elif remote == 'yes':
        path2HeRuCellList = '../Example_cells' 
        path2Base = '../Base/'
    prog = 'minixy_e'
    deltaZ = a_z / 2
    lmp_time = 50
    nprocessor = 2
    if remote == 'no':
        lmp_prog = '/home/heruhe/Downloads/lammps/build/lmp'
    elif remote == 'yes':
        lmp_prog = '/proj/jinxinjx/lammps-2Aug2023/src/lmp_mpi'
    else:
        raise ValueError("the remote should be yes or no.")
    nprint = 2
    
    connect(path2HeRuCellList, path2Base, prog, deltaZ, lmp_time, nprocessor, lmp_prog, nprint, a_z, n_x)
    
    t2 = time.time()
    th, tmin, ts = t2hms(t2 - t1)
    print("Total time: %d h %d min %.1f s" % (th, tmin, ts))






