#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  5 18:28:56 2019

@author: jin

Create crystals

Path:../Modules/Crystal/Create.py

Last modified: 2019-11-28
"""

import numpy as np
import time

class Crystal:
    """
    class about create crystals based on the primitive cell
    """
    
    def __init__(self):
        """
        self.lattice: Lattice parameters along x, y and z direction
        self.element: Targer elements
        self.extend: cell numbers along each direction

        self.xyz: Target atoms coordinates matrix
        self.TotNum: Total target atom numbers
        self.boxsize: Target box size
        self.elementNum: The corresponding number for the element (used in the .in file)
        self.elementCountList: The list of total atoms number belonging to certain element
        self.elementIndexList: The list of element index
        """
        
        self.lattice_x, self.lattice_y, self.lattice_z = 2, 2, 2
        self.element = {"H": 1}
        self.elementSequence = []
        self.extend_x, self.extend_y, self.extend_z = 1, 1, 1
        
        self.xyz = np.empty([1,3])
        self.TotNum = 0
        self.boxsize_x, self.boxsize_y, self.boxsize_z = 2, 2, 2
        
        self.elementNum = [] #Not used
        self.elementCountList = [] #invoked in ElementInfo
        self.elementIndexList = [] #invoked in ElementInfo
        
    def __str__(self):
        ToSay = "Lattice parameter: %.3f %.3f %.3f\n"%(self.lattice_x, self.lattice_y, self.lattice_z)
        
        ToSay = ToSay + 'Elements:'
        for i in self.element:
            ToSay = ToSay + " %s"%i
        ToSay = ToSay + '\n'
        ToSay = ToSay + "Cell numbers: %d %d %d\n"%(self.extend_x, self.extend_y, self.extend_z)
        
        return(ToSay)
        
    def CellPrimitive(self, locin):
        """
        Create the cell by reading a file
        containing the basic information of the target atom location
        """
        
        coord=[]
        elementSequence=[]
        
        with open(locin, 'r') as op:
            #Read the first line by op.readline()
            #Skip the second line by next(op)
            #Then, read the element symbol and the coordinates ratio
            self.TotNum = int(op.readline())
            next(op)
            for line in op:
                linesplit = line.split(' ')
                if linesplit[0] in self.element:
                    elementSequence.append(linesplit[0])
                    coord.append([])
                    coord[-1].append(float(linesplit[1]))
                    coord[-1].append(float(linesplit[2]))
                    coord[-1].append(float(linesplit[3]))        
    
        self.elementSequence = np.array(elementSequence)
        #Calculate the actual coordinates
        #by multiply the coordinate ratio with the lattice parameter
        self.xyz = np.array(coord)
        self.xyz[:,0] *= self.lattice_x
        self.xyz[:,1] *= self.lattice_y
        self.xyz[:,2] *= self.lattice_z
        
        self.boxsize_x, self.boxsize_y, self.boxsize_z = self.lattice_x, self.lattice_y, self.lattice_z
        
    def Extend(self):
        """
        Extend the call according to self.extend
        Note:
            This method is much quicker than concatenate cells.
        """
        
        if (self.extend_x < 1 or self.extend_y < 1 or self.extend_z < 1):
            raise ValueError("The cell cannot be shrinked, check the extend value!")
        elif (not isinstance(self.extend_x, int) or not isinstance(self.extend_y, int) or not isinstance(self.extend_z, int)):
            raise ValueError("The extend number must be the integer type!")
        else:
            if self.extend_x > 1:
                #Note:
                #If use "xyzT_o1=self.xyz" instead of "xyzT_o1=self.xyz * 1"
                #The code "xyzT_o1[:, 0] += self.boxsize_x*jx" will directly change self.xyz
                xyzT_o1=self.xyz * 1
                tmp_cell = np.empty([len(self.xyz) * self.extend_x, 3])
                tmp_cell[0: len(self.xyz)] = self.xyz
                tmp_sequence = []
                tmp_sequence.append(self.elementSequence)
                for jx in range(self.extend_x):
                    if jx>0:
                        xyzT_o1[:, 0] += self.boxsize_x
                        tmp_cell[jx * len(self.xyz) : (jx + 1) * len(self.xyz)] = xyzT_o1
                        tmp_sequence.append(self.elementSequence)
                        
                self.xyz = tmp_cell
                tmp_sequence = np.array(tmp_sequence)
                self.elementSequence = tmp_sequence
                self.elementSequence = self.elementSequence.flatten()

                        
            if self.extend_y > 1:
                xyzT_o2 = self.xyz * 1
                tmp_cell = np.empty([len(self.xyz) * self.extend_y, 3])
                tmp_cell[0: len(self.xyz)] = self.xyz
                tmp_sequence = []
                tmp_sequence.append(self.elementSequence)
                for jy in range(self.extend_y):
                    if jy > 0 :
                        xyzT_o2[:, 1] += self.boxsize_y
                        tmp_cell[jy * len(self.xyz) : (jy + 1) * len(self.xyz)] = xyzT_o2
                        tmp_sequence.append(self.elementSequence)
                        
                self.xyz = tmp_cell
                tmp_sequence = np.array(tmp_sequence)
                self.elementSequence = tmp_sequence
                self.elementSequence = self.elementSequence.flatten()

            if self.extend_z > 1:
                xyzT_o3 = self.xyz * 1
                tmp_cell = np.empty([len(self.xyz) * self.extend_z, 3])
                tmp_cell[0: len(self.xyz)] = self.xyz
                tmp_sequence = []
                tmp_sequence.append(self.elementSequence)
                for jz in range(self.extend_z):
                    if jz > 0 :
                        xyzT_o3[:, 2] += self.boxsize_z
                        tmp_cell[jz * len(self.xyz) : (jz + 1) * len(self.xyz)] = xyzT_o3
                        tmp_sequence.append(self.elementSequence)
                        
                self.xyz = tmp_cell
                tmp_sequence = np.array(tmp_sequence)
                self.elementSequence = tmp_sequence
                self.elementSequence = self.elementSequence.flatten()
                        
            self.boxsize_x *= self.extend_x
            self.boxsize_y *= self.extend_y
            self.boxsize_z *= self.extend_z
            self.TotNum *= self.extend_x * self.extend_y * self.extend_z
            
    def ElementInfo(self):
        """
        Count the element number (self.elementCountList)
        and the element index (self.elementIndexList) in self.xyz.
        Slave of RDA method and remove_atoms.
        """
        self.elementCountList = [] #invoked in ElementInfo
        self.elementIndexList = [] #invoked in ElementInfo
        
        for i in self.element:
            self.elementIndexList.append([])
                    
        for i in range(len(self.elementSequence)):
            self.elementIndexList[int(self.element[self.elementSequence[i]] - 1)].append(i)
        
        for i in self.element:
            self.elementIndexList[int(self.element[i] - 1)] = np.array(self.elementIndexList[int(self.element[i] - 1)])
            self.elementCountList.append(len(self.elementIndexList[int(self.element[i] - 1)]))
    
    def remove_atoms(self, atom_percent, seed = 1):
        """
        Remove atoms according to the percentage given in atom_percent.
        (different elements can have different percentage)

        Parameters
        ----------
        atom_percent : list of 
                       percentage of atoms will be kept, others will be removed
        seed : The random seed

        Returns
        -------
        None.

        """
        
        #Necessary to update number of atoms in each element and corresponding index
        self.ElementInfo()
        
        np.random.seed(seed)
        xyz_tmp = []
        tmp_sequence = [] 
        
        atom_number = [] #Remanining number of atoms after remove
        atom_index = [] #Remaining index of atoms after remove
        
        #Determine the number of atoms and index of atoms after remove (according to element order)
        for i in range(len(self.elementCountList)):
            atom_left = int(self.elementCountList[i]*atom_percent[i]+0.5)
            
            if atom_left > self.elementCountList[i]:
                raise ValueError('More atoms in the target after remove, something wrong!')
                
            atom_number.append(atom_left)
            
            np.random.shuffle(self.elementIndexList[i])
            
            for j in range(atom_left):
                atom_index.append(self.elementIndexList[i][j])
        
        atom_index = np.array(atom_index)
        
    
        #Create new list of atomic coordinates
        for i in range(len(atom_index)):
                xyz_tmp.append(self.xyz[atom_index[i]])
                tmp_sequence.append(self.elementSequence[atom_index[i]])
        
        xyz_tmp = np.array(xyz_tmp)
        tmp_sequence = np.array(tmp_sequence)
        
        self.xyz = xyz_tmp
        self.elementSequence = tmp_sequence
        
        self.TotNum = len(self.xyz)
                    
    def RDA_A(self, rda, seed = 1, CompPrecision = 1e3):
        """
        Introduce RDA into the cell. (Method A)
        rda: The percentage of randomly displaced atoms;
        seed: The random shuffle seed;
        ComPrecision: The precision used to compare float numbers,
                      e.g., 1e3 -> 0.0013 == 0.0009
        Note:
            This method displaces all elements.
            So if one wants to check the displacement of one element,
            one may find that the element target is not located in some region,
            since other element targets occupied those regions
        """
        
        """Step-1: Create a 2D point grid"""
        #The point number along x and y direction is int(np.sqrt(self.TotNum)) + 1
        #Grid_x = [0,1,2,3,...]
        #Grid_y = [0,1,2,3,...]
        #Grid_xy = [[0,0],[0,1],[0,2],[1,0],[1,1],[1,2],[2,0],[2,1],[2,2]]
        #The new 2D location of randomly displaced atom will be chosen from Grid_xy 
        GridPoint = int(np.sqrt(self.TotNum))
        Grid_x = []
        Grid_y = []
        Grid_xy = []
        
        Grid_x = [ i * (self.boxsize_x/GridPoint) for i in range(GridPoint + 1)]
        Grid_y = [ i * (self.boxsize_y/GridPoint) for i in range(GridPoint + 1)]

        Grid_x = np.array(Grid_x)
        Grid_y = np.array(Grid_y)

        for i in range(len(Grid_x)):
            for j in range(len(Grid_y)):
                Grid_xy.append([Grid_x[i], Grid_y[j]])
            
        Grid_xy = np.array(Grid_xy)
        print("Initial grid points: ")
        print(Grid_xy)
        print("Length of grid points:",len(Grid_xy))
        
        #np.random.seed: random number seed;
        #CompPrecision: float number comparison precision, e.g., if it is 1e0, then 2.1 == 1.9;
        """Step-2: Select the RDA list (RandNum) by shuffle the total atom list (TotNumList);"""
        #Also select the list of non-displaced atoms (NonRandNumList);
        #Obtain the X-Y coordinates of the non-displaced atoms (NonRandCoord2D)
        np.random.seed(seed)
        TotNumList = np.arange(0, self.TotNum, 1)
        RandNum = int(self.TotNum * rda / 100 + 0.5)
        np.random.shuffle(TotNumList)
        RandNumList = TotNumList[ : RandNum]
        RandNumList = np.sort(RandNumList)
        NonRandNumList = np.setdiff1d(TotNumList, RandNumList)
        print("RDA:", rda, "Displaced number:", RandNum)
        
        if rda < 100:
            NonRandCoord2D = []
            for i in range(len(NonRandNumList)):
                NonRandCoord2D.append(self.xyz[i, 0 : 2])
            NonRandCoord2D = np.array(NonRandCoord2D)
        
        #Compare the coordinates in NonRandCoord2D and in Grid_xy (represented by Grid_x and Grid_y);
        #If the value in Grid_xy is found in NonRandCoord2D,
        #then delete that value.
        """Step-3: In order to obtain a new Grid_xy with no point occupied by the non-displaced atoms"""
        #Index(Grid_xy) = Index(Grid_x) * len(Grid_x) + Index(Grid_y)
        #Note: 
        #   If the value of CompPrecision, then the similar part between Grid_xy and NonRandCoord2D is higher.
        if rda < 100:
            NonRandCoord2D = NonRandCoord2D * CompPrecision
            NonRandCoord2D = NonRandCoord2D + 0.5
            NonRandCoord2D = NonRandCoord2D.astype(int)
            Grid_x = Grid_x * CompPrecision
            Grid_x = Grid_x + 0.5
            Grid_x = Grid_x.astype(int)
            Grid_y = Grid_y * CompPrecision
            Grid_y = Grid_y + 0.5
            Grid_y = Grid_y.astype(int)
        
            PopOutIndex=[]
                
            for i in range(len(Grid_x)):
                if Grid_x[i] in NonRandCoord2D[:, 0]:
                    for j in range(len(Grid_y)):
                        if Grid_y[j] in NonRandCoord2D[:, 1]:
                            PopOutIndex.append(i * len(Grid_x) + j)
         
            Grid_xy = np.delete(Grid_xy, PopOutIndex, 0)
           
        print("New grid points: ")
        print(Grid_xy)
        print("Length of grid points:",len(Grid_xy))
        
        """Step-4: Move the randomly displaced atoms to the location indicated by Grid_xy"""
        #without changing the z coordinate
        j = 0 
        np.random.shuffle(Grid_xy)
        if RandNum <= len(Grid_xy):
            for i in range(RandNum):
                #print('Before moved:',i, self.xyz[i])
                self.xyz[RandNumList[i], 0 : 2] = Grid_xy[j]
                #print('After moved:',i, self.xyz[i])
                j += 1
        else:
            raise ValueError("More displaced atoms than avaliable grid points, Check!")
            
    def RDA_B(self, rda, symbol, Point1D = 0, seed = 1, CompPrecision = 1e3):
        """
        Introduce RDA into the cell. (Method B)
        Only the element selected by the symbol will be displaced.
        rda: The percentage of randomly displaced atoms;
        Point1D: The gird point along one direction,
                 If 0, use a default value, if not 0, use user input value;
        symbol: The symbol of the element to be displaced
        seed: The random shuffle seed;
        ComPrecision: The precision used to compare float numbers,
                      e.g., 1e3 -> 0.0013 == 0.0009
        """
        
        symbolValue = int(self.element[symbol] - 1)
        """Step-1: Create a 2D point grid"""
        #The point number along x and y direction is int(np.sqrt(self.TotNum)) + 1
        #Grid_x = [0,1,2,3,...]
        #Grid_y = [0,1,2,3,...]
        #Grid_xy = [[0,0],[0,1],[0,2],[1,0],[1,1],[1,2],[2,0],[2,1],[2,2]]
        #The new 2D location of randomly displaced atom will be chosen from Grid_xy 
        if Point1D == 0:
            GridPoint = int(np.sqrt(self.elementCountList[symbolValue]))
        else:
            GridPoint = Point1D
            
        Grid_x = []
        Grid_y = []
        Grid_xy = []
        
        Grid_x = [ i * (self.boxsize_x/GridPoint) for i in range(GridPoint + 1)]
        Grid_y = [ i * (self.boxsize_y/GridPoint) for i in range(GridPoint + 1)]

        Grid_x = np.array(Grid_x)
        Grid_y = np.array(Grid_y)

        for i in range(len(Grid_x)):
            for j in range(len(Grid_y)):
                Grid_xy.append([Grid_x[i], Grid_y[j]])
            
        Grid_xy = np.array(Grid_xy)
#        print("Initial grid points: ")
#        print(Grid_xy)
#        print("Length of grid points:",len(Grid_xy))

        #np.random.seed: random number seed;
        #CompPrecision: float number comparison precision, e.g., if it is 1e0, then 2.1 == 1.9;
        """Step-2: Select the RDA list (RandNum) by shuffle the total atom list (TotNumList);"""
        #Also select the list of non-displaced atoms (NonRandNumList);
        #Obtain the X-Y coordinates of the non-displaced atoms (NonRandCoord2D)
        np.random.seed(seed)
        TotNumList = np.arange(0, self.elementCountList[symbolValue], 1)
        RandNum = int(self.elementCountList[symbolValue] * rda / 100 + 0.5)
        np.random.shuffle(TotNumList)
        RandNumList = TotNumList[ : RandNum]
        RandNumList = np.sort(RandNumList)
        NonRandNumList = np.setdiff1d(TotNumList, RandNumList)
        
        if rda < 100:
            NonRandCoord2D = []
            for i in range(len(NonRandNumList)):
                NonRandCoord2D.append(self.xyz[i, 0 : 2])
            NonRandCoord2D = np.array(NonRandCoord2D)
            
        #Compare the coordinates in NonRandCoord2D and in Grid_xy (represented by Grid_x and Grid_y);
        #If the value in Grid_xy is found in NonRandCoord2D,
        #then delete that value.
        """Step-3: In order to obtain a new Grid_xy with no point occupied by the non-displaced atoms"""
        #Index(Grid_xy) = Index(Grid_x) * len(Grid_x) + Index(Grid_y)
        #Note: 
        #   If the value of CompPrecision, then the similar part between Grid_xy and NonRandCoord2D is higher.
        if rda < 100:
            NonRandCoord2D = NonRandCoord2D * CompPrecision
            NonRandCoord2D = NonRandCoord2D + 0.5
            NonRandCoord2D = NonRandCoord2D.astype(int)
            Grid_x = Grid_x * CompPrecision
            Grid_x = Grid_x + 0.5
            Grid_x = Grid_x.astype(int)
            Grid_y = Grid_y * CompPrecision
            Grid_y = Grid_y + 0.5
            Grid_y = Grid_y.astype(int)
        
            PopOutIndex=[]
                
            for i in range(len(Grid_x)):
                if Grid_x[i] in NonRandCoord2D[:, 0]:
                    for j in range(len(Grid_y)):
                        if Grid_y[j] in NonRandCoord2D[:, 1]:
                            PopOutIndex.append(i * len(Grid_x) + j)
         
            Grid_xy = np.delete(Grid_xy, PopOutIndex, 0)
           
#        print("New grid points: ")
#        print(Grid_xy)
#        print("Length of grid points:",len(Grid_xy))
        
        """Step-4: Move the randomly displaced atoms to the location indicated by Grid_xy"""
        #without changing the z coordinate
        j = 0 
        np.random.shuffle(Grid_xy)
        if RandNum <= len(Grid_xy):
            for i in range(RandNum):
                #print('Before moved:',i, self.xyz[i])
                self.xyz[self.elementIndexList[symbolValue][RandNumList[i]], 0 : 2] = Grid_xy[j]
                #print('After moved:',i, self.xyz[i])
                j += 1
        else:
            raise ValueError("More displaced atoms than avaliable grid points, Check!")

    def RDA_C(self, rda, symbol, dist, z_flag = 0, seed = 1):
        """
        Introduce RDA into the cell. (Method C)
        Only the element selected by the symbol will be displaced.
        rda: The percentage of randomly displaced atoms;
        symbol: The symbol of the element to be displaced;
        dist: The range of atom to be randomly displaced along each dimension
        z_flag: If z_flag = 0, the atom is allowed to be displaced along z direction
                If z_flat = 1, tha atom is not allowed to move along z direction;
        seed: The random shuffle seed;
        """
        
        symbolValue = int(self.element[symbol] - 1)
        
        """Step-1: Select the RDA list (RandNum) by shuffle the total atom list (TotNumList);"""
        np.random.seed(seed)
        TotNumList = np.arange(0, self.elementCountList[symbolValue], 1)
        RandNum = int(self.elementCountList[symbolValue] * rda / 100 + 0.5)
        np.random.shuffle(TotNumList)
        RandNumList = TotNumList[ : RandNum]
        RandNumList = np.sort(RandNumList)
        
        """Step-2: Creat list of randomly displaced distance"""
        #from -1*(dist)/2 to (dist)/2
        xR=np.empty(RandNum)
        yR=np.empty(RandNum)
        
        for i in range(RandNum):
            xR[i]=(dist*(np.random.random_sample()-0.5))
            yR[i]=(dist*(np.random.random_sample()-0.5))
        
        if z_flag == 0:
            zR=np.empty(RandNum)
            for i in range(RandNum):
                zR[i]=(dist*(np.random.random_sample()-0.5))
                
        """Step-3: Move the randomly displace atoms"""
        for i in range(RandNum):
            self.xyz[self.elementIndexList[symbolValue][RandNumList[i]], 0] += xR[i]
            if self.xyz[self.elementIndexList[symbolValue][RandNumList[i]], 0] < 0:
                self.xyz[self.elementIndexList[symbolValue][RandNumList[i]], 0] += self.boxsize_x
            if self.xyz[self.elementIndexList[symbolValue][RandNumList[i]], 0] >=  self.boxsize_x:
                self.xyz[self.elementIndexList[symbolValue][RandNumList[i]], 0] -= self.boxsize_x
                
            self.xyz[self.elementIndexList[symbolValue][RandNumList[i]], 1] += yR[i]
            if self.xyz[self.elementIndexList[symbolValue][RandNumList[i]], 1] < 0:
                self.xyz[self.elementIndexList[symbolValue][RandNumList[i]], 1] += self.boxsize_y
            if self.xyz[self.elementIndexList[symbolValue][RandNumList[i]], 1] >=  self.boxsize_y:
                self.xyz[self.elementIndexList[symbolValue][RandNumList[i]], 1] -= self.boxsize_y

            if z_flag == 0:
                self.xyz[self.elementIndexList[symbolValue][RandNumList[i]], 2] += zR[i]
                if self.xyz[self.elementIndexList[symbolValue][RandNumList[i]], 2] < 0:
                    self.xyz[self.elementIndexList[symbolValue][RandNumList[i]], 2] += self.boxsize_z
                if self.xyz[self.elementIndexList[symbolValue][RandNumList[i]], 2] >=  self.boxsize_z:
                    self.xyz[self.elementIndexList[symbolValue][RandNumList[i]], 2] -= self.boxsize_z

    def RDA_C2(self, rda, dist, z_flag = 0, seed = 1):
        """
        Introduce RDA into the cell. (Method C)
        Displace all type of atoms
        rda: The percentage of randomly displaced atoms;
        dist: The range of atom to be randomly displaced along each dimension
        z_flag: If z_flag = 0, the atom is allowed to be displaced along z direction
                If z_flat = 1, tha atom is not allowed to move along z direction;
        seed: The random shuffle seed;
        """
        
        """Step-1: Select the RDA list (RandNum) by shuffle the total atom list (TotNumList);"""
        np.random.seed(seed)
        TotNumList = np.arange(0, self.TotNum, 1)
        RandNum = int(self.TotNum * rda / 100 + 0.5)
        np.random.shuffle(TotNumList)
        RandNumList = TotNumList[ : RandNum]
        RandNumList = np.sort(RandNumList)
        
        """Step-2: Creat list of randomly displaced distance"""
        #from -1*(dist)/2 to (dist)/2
        xR=np.empty(RandNum)
        yR=np.empty(RandNum)
        
        for i in range(RandNum):
            xR[i]=(dist*(np.random.random_sample()-0.5))
            yR[i]=(dist*(np.random.random_sample()-0.5))
        
        if z_flag == 0:
            zR=np.empty(RandNum)
            for i in range(RandNum):
                zR[i]=(dist*(np.random.random_sample()-0.5))
                
        """Step-3: Move the randomly displaced atoms"""
        for i in range(RandNum):
            self.xyz[RandNumList[i], 0] += xR[i]
            if self.xyz[RandNumList[i], 0] < 0:
                self.xyz[RandNumList[i], 0] += self.boxsize_x
            if self.xyz[RandNumList[i], 0] >=  self.boxsize_x:
                self.xyz[RandNumList[i], 0] -= self.boxsize_x
                
            self.xyz[RandNumList[i], 1] += yR[i]
            if self.xyz[RandNumList[i], 1] < 0:
                self.xyz[RandNumList[i], 1] += self.boxsize_y
            if self.xyz[RandNumList[i], 1] >=  self.boxsize_y:
                self.xyz[RandNumList[i], 1] -= self.boxsize_y

            if z_flag == 0:
                self.xyz[RandNumList[i], 2] += zR[i]
                if self.xyz[RandNumList[i], 2] < 0:
                    self.xyz[RandNumList[i], 2] += self.boxsize_z
                if self.xyz[RandNumList[i], 2] >=  self.boxsize_z:
                    self.xyz[RandNumList[i], 2] -= self.boxsize_z
                    
    def RDA_C3(self, rda, distx, disty, distz, z_flag = 0, seed = 1):
        """
        Similar to function RDA_C2,
        allow to have different displacement range in x, y and z directions:
            distx, disty and distz
        """
        
        """Step-1: Select the RDA list (RandNum) by shuffle the total atom list (TotNumList);"""
        np.random.seed(seed)
        TotNumList = np.arange(0, self.TotNum, 1)
        RandNum = int(self.TotNum * rda / 100 + 0.5)
        np.random.shuffle(TotNumList)
        RandNumList = TotNumList[ : RandNum]
        RandNumList = np.sort(RandNumList)
        
        """Step-2: Creat list of randomly displaced distance"""
        #from -1*(dist)/2 to (dist)/2
        xR=np.empty(RandNum)
        yR=np.empty(RandNum)
        
        for i in range(RandNum):
            xR[i]=(distx*(np.random.random_sample()-0.5))
            yR[i]=(disty*(np.random.random_sample()-0.5))
        
        if z_flag == 0:
            zR=np.empty(RandNum)
            for i in range(RandNum):
                zR[i]=(distz*(np.random.random_sample()-0.5))
                
        """Step-3: Move the randomly displaced atoms"""
        for i in range(RandNum):
            self.xyz[RandNumList[i], 0] += xR[i]
            if self.xyz[RandNumList[i], 0] < 0:
                self.xyz[RandNumList[i], 0] += self.boxsize_x
            if self.xyz[RandNumList[i], 0] >=  self.boxsize_x:
                self.xyz[RandNumList[i], 0] -= self.boxsize_x
                
            self.xyz[RandNumList[i], 1] += yR[i]
            if self.xyz[RandNumList[i], 1] < 0:
                self.xyz[RandNumList[i], 1] += self.boxsize_y
            if self.xyz[RandNumList[i], 1] >=  self.boxsize_y:
                self.xyz[RandNumList[i], 1] -= self.boxsize_y

            if z_flag == 0:
                self.xyz[RandNumList[i], 2] += zR[i]
                if self.xyz[RandNumList[i], 2] < 0:
                    self.xyz[RandNumList[i], 2] += self.boxsize_z
                if self.xyz[RandNumList[i], 2] >=  self.boxsize_z:
                    self.xyz[RandNumList[i], 2] -= self.boxsize_z    
                    
    def RDA_C4(self, rda, x, y, z, z_flag = 0, seed = 1):
        """
        Similar to function RDA_C3,
        The atoms will be displaced according to a fixed distance:
            x, y and z
        """
        
        """Step-1: Select the RDA list (RandNum) by shuffle the total atom list (TotNumList);"""
        np.random.seed(seed)
        TotNumList = np.arange(0, self.TotNum, 1)
        RandNum = int(self.TotNum * rda / 100 + 0.5)
        np.random.shuffle(TotNumList)
        RandNumList = TotNumList[ : RandNum]
        RandNumList = np.sort(RandNumList)
        
        """Step-2: Creat list of randomly displaced distance"""
        #from -1*(dist)/2 to (dist)/2
        xR = x
        yR = y
        zR = z
                
        """Step-3: Move the randomly displaced atoms"""
        for i in range(RandNum):
            self.xyz[RandNumList[i], 0] += xR
            if self.xyz[RandNumList[i], 0] < 0:
                self.xyz[RandNumList[i], 0] += self.boxsize_x
            if self.xyz[RandNumList[i], 0] >=  self.boxsize_x:
                self.xyz[RandNumList[i], 0] -= self.boxsize_x
                
            self.xyz[RandNumList[i], 1] += yR
            if self.xyz[RandNumList[i], 1] < 0:
                self.xyz[RandNumList[i], 1] += self.boxsize_y
            if self.xyz[RandNumList[i], 1] >=  self.boxsize_y:
                self.xyz[RandNumList[i], 1] -= self.boxsize_y

            if z_flag == 0:
                self.xyz[RandNumList[i], 2] += zR
                if self.xyz[RandNumList[i], 2] < 0:
                    self.xyz[RandNumList[i], 2] += self.boxsize_z
                if self.xyz[RandNumList[i], 2] >=  self.boxsize_z:
                    self.xyz[RandNumList[i], 2] -= self.boxsize_z  
                    
    def RDA_McChasy(self, rda_file, dist, z_flag = 0, seed = 1):
        """
        Introduce the RDA depth profile like the one used in McChasy.
        rda_file: 
            Path to the RDA depth profile;
            The first row is the depth [nm];
            The second row is the RDA [%];
        dist: The range of atom to be randomly displaced along each dimension
        z_flag: If z_flag = 0, the atom is allowed to be displaced along z direction
                If z_flat = 1, tha atom is not allowed to move along z direction;
        seed: The random shuffle seed;
        """
        
        """Step-1: Read the RDA depth profile"""
        Depth_RDA = []
        RDA = []
        
        with open(rda_file, 'r') as op:
            next(op)
            for line in op:
                linesplit = line.split()
                Depth_RDA.append(float(linesplit[0])*10)
                RDA.append(float(linesplit[-1]))
           
        if Depth_RDA[-1] >= self.boxsize_z:
            raise ValueError("The maximum depth (%.3f) in the RDA profile should not be larger than the target depth (%.3f)!"\
                             %(Depth_RDA[-1], self.boxsize_z))
        
        Depth_RDA.append(self.boxsize_z)
        Depth_RDA, RDA = np.array(Depth_RDA), np.array(RDA)
        
        """Step-2: Select atoms in the corresponding depth region"""
        #Seperate the atom into different groups
        #according to the depth
        AtomIndexList = []
        for i in range(len(RDA)):
            AtomIndexList.append([])
            
        for i in range(len(self.xyz)):            
            for j in range(1, len(Depth_RDA)):
                if self.xyz[i,2] >= Depth_RDA[j-1] and self.xyz[i,2] <= Depth_RDA[j]:
                    AtomIndexList[j-1].append(i)
        
        for i in range(len(AtomIndexList)):
            AtomIndexList[i] = np.array(AtomIndexList[i])
        
        """Step-3: Displace atoms according to the RDA profile"""
        np.random.seed(seed)
        
        for k in range(len(RDA)):
            #Obtain the index of atoms which needs to be displaced.
            RandNum = int(len(AtomIndexList[k]) * RDA[k] / 100 + 0.5)
            np.random.shuffle(AtomIndexList[k])
            RandNumList = AtomIndexList[k][ : RandNum]
            RandNumList = np.sort(RandNumList)
            
            xR=np.empty(RandNum)
            yR=np.empty(RandNum)
        
            #Prepare the displacement distance
            for i in range(RandNum):
                xR[i]=(dist*(np.random.random_sample()-0.5))
                yR[i]=(dist*(np.random.random_sample()-0.5))
        
            if z_flag == 0:
                zR=np.empty(RandNum)
                for i in range(RandNum):
                    zR[i]=(dist*(np.random.random_sample()-0.5))
                    
            for i in range(RandNum):
                self.xyz[RandNumList[i], 0] += xR[i]
                if self.xyz[RandNumList[i], 0] < 0:
                    self.xyz[RandNumList[i], 0] += self.boxsize_x
                if self.xyz[RandNumList[i], 0] >=  self.boxsize_x:
                    self.xyz[RandNumList[i], 0] -= self.boxsize_x
                
                self.xyz[RandNumList[i], 1] += yR[i]
                if self.xyz[RandNumList[i], 1] < 0:
                    self.xyz[RandNumList[i], 1] += self.boxsize_y
                if self.xyz[RandNumList[i], 1] >=  self.boxsize_y:
                    self.xyz[RandNumList[i], 1] -= self.boxsize_y
                
                #The displacement along the z direction is different with
                #the traditional RDA method.
                #Here, the Depth_RDA is used instead of self.boxsize_z
                if z_flag == 0:
                    self.xyz[RandNumList[i], 2] += zR[i]
                if self.xyz[RandNumList[i], 2] < Depth_RDA[k]:
                    self.xyz[RandNumList[i], 2] += (Depth_RDA[k+1] - Depth_RDA[k])
                if self.xyz[RandNumList[i], 2] >=  Depth_RDA[k+1]:
                    self.xyz[RandNumList[i], 2] -= (Depth_RDA[k+1] - Depth_RDA[k])
                    
    def RDA_McChasy2(self, rda_file, distx, disty, distz, z_flag = 0, seed = 1):
        """
        Similar to RDA_McChasy().
        allow to have different displacement range in x, y and z directions:
            distx, disty and distz
        """
    
        """Step-1: Read the RDA depth profile"""
        Depth_RDA = []
        RDA = []
        
        with open(rda_file, 'r') as op:
            next(op)
            for line in op:
                linesplit = line.split()
                Depth_RDA.append(float(linesplit[0])*10)
                RDA.append(float(linesplit[-1]))
           
        if Depth_RDA[-1] >= self.boxsize_z:
            raise ValueError("The maximum depth (%.3f) in the RDA profile should not be larger than the target depth (%.3f)!"\
                             %(Depth_RDA[-1], self.boxsize_z))
        
        Depth_RDA.append(self.boxsize_z)
        Depth_RDA, RDA = np.array(Depth_RDA), np.array(RDA)
        
        """Step-2: Select atoms in the corresponding depth region"""
        #Seperate the atom into different groups
        #according to the depth
        AtomIndexList = []
        for i in range(len(RDA)):
            AtomIndexList.append([])
            
        for i in range(len(self.xyz)):            
            for j in range(1, len(Depth_RDA)):
                if self.xyz[i,2] >= Depth_RDA[j-1] and self.xyz[i,2] <= Depth_RDA[j]:
                    AtomIndexList[j-1].append(i)
        
        for i in range(len(AtomIndexList)):
            AtomIndexList[i] = np.array(AtomIndexList[i])
        
        """Step-3: Displace atoms according to the RDA profile"""
        np.random.seed(seed)
        
        for k in range(len(RDA)):
            #Obtain the index of atoms which needs to be displaced.
            RandNum = int(len(AtomIndexList[k]) * RDA[k] / 100 + 0.5)
            np.random.shuffle(AtomIndexList[k])
            RandNumList = AtomIndexList[k][ : RandNum]
            RandNumList = np.sort(RandNumList)
            
            xR=np.empty(RandNum)
            yR=np.empty(RandNum)
        
            #Prepare the displacement distance
            for i in range(RandNum):
                xR[i]=(distx*(np.random.random_sample()-0.5))
                yR[i]=(disty*(np.random.random_sample()-0.5))
        
            if z_flag == 0:
                zR=np.empty(RandNum)
                for i in range(RandNum):
                    zR[i]=(distz*(np.random.random_sample()-0.5))
                    
            for i in range(RandNum):
                self.xyz[RandNumList[i], 0] += xR[i]
                if self.xyz[RandNumList[i], 0] < 0:
                    self.xyz[RandNumList[i], 0] += self.boxsize_x
                if self.xyz[RandNumList[i], 0] >=  self.boxsize_x:
                    self.xyz[RandNumList[i], 0] -= self.boxsize_x
                
                self.xyz[RandNumList[i], 1] += yR[i]
                if self.xyz[RandNumList[i], 1] < 0:
                    self.xyz[RandNumList[i], 1] += self.boxsize_y
                if self.xyz[RandNumList[i], 1] >=  self.boxsize_y:
                    self.xyz[RandNumList[i], 1] -= self.boxsize_y
                
                #The displacement along the z direction is different with
                #the traditional RDA method.
                #Here, the Depth_RDA is used instead of self.boxsize_z
                if z_flag == 0:
                    self.xyz[RandNumList[i], 2] += zR[i]
                if self.xyz[RandNumList[i], 2] < Depth_RDA[k]:
                    self.xyz[RandNumList[i], 2] += (Depth_RDA[k+1] - Depth_RDA[k])
                if self.xyz[RandNumList[i], 2] >=  Depth_RDA[k+1]:
                    self.xyz[RandNumList[i], 2] -= (Depth_RDA[k+1] - Depth_RDA[k])

    def Write_in(self, locout):
        """
        Write the cell with the .in format
        """
        
#        for i in range(len(self.elementSequence)):
#            self.elementNum.append(self.element[self.elementSequence[i]])
            
        with open(locout,'w') as wt:
            wt.write('%d\n'%(self.TotNum))
            wt.write('Lattice="%.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f"\n'%(self.boxsize_x,0.0,0.0,0.0,self.boxsize_y,0.0,0.0,0.0,self.boxsize_z))
            for i in range(len(self.xyz)):
                wt.write('%s %.5f %.5f %.5f %s\n'%(self.elementSequence[i],self.xyz[i,0],self.xyz[i,1],self.xyz[i,2],self.element[self.elementSequence[i]]))

    
if __name__ == '__main__':
    print("Hello")
    #start = time.time()
    
    """Test"""
#    UO2A = Crystal()
#    UO2A.lattice_x, UO2A.lattice_y, UO2A.lattice_z = 5.471, 5.471, 5.471
#    UO2A.element = {"U": 1, "O": 2}
#    UO2A.extend_x, UO2A.extend_y, UO2A.extend_z = 40, 40, 40
#    print(UO2A)
#    
#    PrimitCellLoc = '/home/jin/Licument/RBS_Simu/CellCell/PrimUO2.xyz'
#    UO2A.CellPrimitive(PrimitCellLoc)
#    
#    UO2A.Extend()
#    UO2A.ElementInfo()
#    time_rda = time.time()
##    UO2A.RDA(100)
#    #UO2A.RDA_B(100, 'U', Point1D = 2962)
#    UO2A.RDA_C(100, 'U', 5.471, z_flag = 1)
#
##    
#    end1 = time.time()
#    TestOut = 'UO2A_test.in'
#    UO2A.Write_in(TestOut)
#    
#    end2 = time.time()
#    print("Time spent on RDA:", end1 - time_rda)
#    print("Time spent on write:", end2 - end1)
#    print("Total time elapsed:", end2 - start)
    """Test"""
    
    """For si 2019/09/27 A"""
#    Si_A = Crystal()
#    Si_A.lattice_x, Si_A.lattice_y, Si_A.lattice_z = 5.4307, 5.4307, 5.4307
#    Si_A.element = {"Si": 1}
#    Si_A.extend_x, Si_A.extend_y, Si_A.extend_z = 1, 1, 1
#    print(Si_A)
#    
#    Si_A_PrimitCellLoc = '/home/jin/Licument/RBS_Simu/CellCell/PrimSi_A.xyz'
#    Si_A.CellPrimitive(Si_A_PrimitCellLoc)
#    Si_A.Extend()
#    Si_A_Out = '/home/jin/Licument/RBS_Simu/Crystal_structure/Si_A/coords1.in'
#    Si_A.Write_in(Si_A_Out)
    """For si 2019/09/27 A"""
    
    """For si 2019/09/27 B"""
#    Si_B = Crystal()
#    Si_B.lattice_x, Si_B.lattice_y, Si_B.lattice_z = 5.4307, 5.4307, 5.4307
#    Si_B.element = {"Si": 1}
#    Si_B.extend_x, Si_B.extend_y, Si_B.extend_z = 1, 1, 1
#    print(Si_B)
#    
#    Si_B_PrimitCellLoc = '/home/jin/Licument/RBS_Simu/Crystal_structure/Si_B/PrimSi_B.xyz'
#    Si_B.CellPrimitive(Si_B_PrimitCellLoc)
#    Si_B.Extend()
#    Si_B_Out = '/home/jin/Licument/RBS_Simu/Crystal_structure/Si_B/coords1.in'
#    Si_B.Write_in(Si_B_Out)
    """For si 2019/09/27 B"""
    
    """Test UO2 2019/10/03"""
#    UO2A = Crystal()
#    UO2A.lattice_x, UO2A.lattice_y, UO2A.lattice_z = 5.471, 5.471, 5.471
#    UO2A.element = {"U": 1, "O": 2}
#    UO2A.extend_x, UO2A.extend_y, UO2A.extend_z = 40, 40, 40
#    print(UO2A)
#    
#    PrimitCellLoc = '/home/jin/Licument/RBS_Simu/CellCell/PrimUO2.xyz'
#    UO2A.CellPrimitive(PrimitCellLoc)
#    
#    UO2A.Extend()
#    #UO2A.ElementInfo()
#    #time_rda = time.time()
#    #UO2A.RDA_C2(100, 5.471, z_flag = 0, seed = 1)
#    #end1 = time.time()
#    UO2A.RDA_McChasy('/home/jin/Licument/Modules/Crystal/RDA_Profile_tmpA.txt', 5.471, z_flag = 0, seed = 1)
#    
#    TestOut = 'UO2A_test.in'
#    UO2A.Write_in(TestOut)
    
#    end2 = time.time()
#    print("Time spent on RDA:", end1 - time_rda)
#    print("Time spent on write:", end2 - end1)
#    print("Total time elapsed:", end2 - start)
    """Test UO2 2019/10/03"""
    
    """For GaAs 2019/10/21"""
#    GaAs_A = Crystal()
#    GaAs_A.lattice_x, GaAs_A.lattice_y, GaAs_A.lattice_z = 5.6532, 5.6532, 5.6532
#    GaAs_A.element = {"Ga": 1, "As": 2}
#    GaAs_A.extend_x, GaAs_A.extend_y, GaAs_A.extend_z = 1, 1, 1
#    print(GaAs_A)
#    
#    GaAs_A_PrimitCellLoc = '/home/jin/Licument/RBS_Simu/Crystal_structure/GaAs_A/PrimGaAs_A.xyz'
#    GaAs_A.CellPrimitive(GaAs_A_PrimitCellLoc)
#    GaAs_A.Extend()
#    GaAs_A_Out = '/home/jin/Licument/RBS_Simu/Crystal_structure/GaAs_A/coords1.in'
#    GaAs_A.Write_in(GaAs_A_Out)
    """For GaAs 2019/10/21"""
    
    """For Fe 2019/11/26 A"""
#    Fe_A = Crystal()
#    Fe_A.lattice_x, Fe_A.lattice_y, Fe_A.lattice_z = 2.8665, 2.8665, 2.8665
#    Fe_A.element = {"Fe": 1}
#    Fe_A.extend_x, Fe_A.extend_y, Fe_A.extend_z = 1, 1, 1
#    print(Fe_A)
#    
#    Fe_A_PrimitCellLoc = '/home/jin/Licument/RBS_Simu/CellCell/PrimFe_A.xyz'
#    Fe_A.CellPrimitive(Fe_A_PrimitCellLoc)
#    Fe_A.Extend()
#    Fe_A_Out = '/home/jin/Licument/RBS_Simu/Crystal_structure/Fe_A/coords1.in'
#    Fe_A.Write_in(Fe_A_Out)
    """For Fe 2019/11/26 A"""

    
