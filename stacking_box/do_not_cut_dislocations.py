#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 25 11:34:02 2022

@author: jinxinjx

    
    __|    
    __|    |
    __|    |    |
    __|    |    |
    __|         |
    __|  _ _ _ _ _ _ _ _ _ _ (Middle of the biggest un-covered region)
    __|
    __|      |
    __|      |
    
    As shown above, find the biggest region that are not covered by any 
    dislocations, then select the middle of that biggest region as the 
    new interface, then move the cell below the interface to top and move
    the cell above the interface to bottom.
"""

from ovito.io import import_file, export_file
from ovito.modifiers import *
from ovito.data import DislocationNetwork

import ovito
import os
import numpy as np
import pandas as pd

def uncut(filein, fileout, binWidth, limit_length):
    """

    Parameters
    ----------
    filein : string
        The path to the input file.
    fileout : string
        The path to the output file
    binWidth : 
        The width of the bin that discritize the cell
    limit_length : 
        The minimum length of the region that is not covered by any dislocations.

    Returns
    -------
    None.
    (adjust the cell boundary)
    """
    
    # Enables activity logging for OVITO’s internal operations
    #ovito.enable_logging()
    
    pipeline = import_file(filein)

    data = pipeline.compute()
    
    """
    1:
      Discretize the cell along the z direction by a list of bin.
      
      We need to provide the width of the bin (may be modified by code),
      and a limit length from the function arguments.
      (If the uncovered region is smaller than this limit_length, we think
       the dislocations are too close and the cell boundaries will not be
       adjusted.)
    """
    
    #Get the lower limit and upper limit of the cell along the z direction
    z_low = data.cell[-1,-1]
    z_high = data.cell[-1,-1] + data.cell[-1,-2]
    z_length = data.cell[-1,-2]
    
    binNumber = int(z_length / binWidth)
    
    #Recalculate the bin width according to the bin number
    binWidth = z_length / binNumber
    
    binList = np.arange(z_low, z_high, binWidth)
    
    """
    2:
        Extract dislocations by using OVITO;
        Find dislocation segments.
        
        If no dislocation is cut by the z plane, flag=0, no operation;
        
        Otherwise, flag = 1, start the operation:
            split the dislocation cut by the z plane;
            move the newly generated dislocations according to the rule
            of periodic bondary condition.
            
        Note:
            if a dislocation is cut by the boundary, some dislocation segments
            calculated by OVITO will be out of the boundary limit.
            (The whole program is based on this property of OVITO)
        
    """
    # Extract dislocation lines from a crystal with BCC structure:
    modifier = DislocationAnalysisModifier()
    modifier.input_crystal_structure = DislocationAnalysisModifier.Lattice.BCC
    pipeline.modifiers.append(modifier)
    
    data = pipeline.compute()
    
    # Print list of dislocation lines:
    # Get the list of the z position of dislocation segments
    # (We only need the minimum and the maximum of the z position)
    flag = 0
    disloc_segment_z=[]
    #print("Found %i dislocation segments" % len(data.dislocations.segments))
    for segment in data.dislocations.segments:
        #print("Segment %i: length=%f, Burgers vector=%s" % (segment.id, segment.length, segment.true_burgers_vector))
        #print(segment.points)
        disz_mini = np.min(segment.points[:,-1])
        disz_max = np.max(segment.points[:,-1])
        
        if disz_mini >= z_low and disz_max <= z_high :
            disloc_segment_z.append([disz_mini, disz_max])
            
        elif disz_mini < z_low and disz_max > z_high :
            disloc_segment_z.append([z_low, z_high])
        
        #Split a dislocation to two, and apply the peridoic boundary condition    
        elif disz_mini < z_low and disz_max <= z_high :
            flag = 1
            disloc_segment_z.append([z_low, disz_max])
            disloc_segment_z.append([z_high - (z_low - disz_mini), z_high])
            
        #Split a dislocation to two, and apply the peridoic boundary condition    
        elif disz_mini >= z_low and disz_max > z_high :
            flag = 1
            disloc_segment_z.append([disz_mini, z_high])
            disloc_segment_z.append([z_low, z_low + (disz_max - z_high)])    
            
        else:
            print("Warning: special condition for disz_mini = %.5f and disz_max = %.5f." % (disz_mini, disz_max))
            print("In file: %s." % filein)
            continue
            
    if flag == 0:
        print("%s: no dislocation is cut by the z planes (no more operation)." % filein)
    else:
        print("%s: Dislocation-cut is spotted." % filein)
        
        """
        3:
            Find the bins that are not covered by dislocation segments along the z direction,
            by looping over the list of dislocation segments (only containing the extreme segments).
            
            The index of the corresponding bins will be stored.
        """
        
        commen_bin_index = []
        for i in range(len(disloc_segment_z)):
            
            commen_bin_index.append([])
            
            for j in range(len(binList)):
                if not (binList[j] >= disloc_segment_z[i][0] and binList[j] <= disloc_segment_z[i][1]):
                    commen_bin_index[-1].append(j)
                    
        for i in range(len(commen_bin_index)):
            commen_bin_index[i] = np.array(commen_bin_index[i])
            
                
        """
        4:
            Compare the list of the bins that are not covered,
            in order to find the commen values.
            
            If no commen value is find, stop the operation;
            Otherwise, continue.
        """
        bin_noCover = commen_bin_index[0]
        
        for i in range(1, len(commen_bin_index)):
            bin_noCover = np.intersect1d(bin_noCover, commen_bin_index[i], assume_unique=True)
            
        
        if len(bin_noCover) < 1:
            print("Warning: all the space is covered by the dislocations! (no more operation)")
        else:
            
            """
            5:
                If the bins are connected with each other,
                they are organized to a chunk.
                
                If the biggest chunk is too small, then no more operation.
            """
            
            bin_chunks = []
            bin_chunks.append([])
            
            #Create the bin chunks
            if len(bin_noCover) == 1:
                bin_chunks[-1].append(bin_noCover[0])
            else:
                for i in range(1, len(bin_noCover)):
                    bin_chunks[-1].append(bin_noCover[i-1])
                    if bin_noCover[i] != (bin_noCover[i-1] + 1):
                        bin_chunks.append([])
                
                bin_chunks[-1].append(bin_noCover[-1])
                            
            
            #Find the biggest chunk
            bin_len_max = 0
            loc_max = 0
            for i in range(len(bin_chunks)):
                if len(bin_chunks[i]) > bin_len_max:
                    bin_len_max = len(bin_chunks[i])
                    loc_max = i
                    
            bin_biggest_chunk = bin_chunks[loc_max]
            
            if bin_len_max * binWidth < limit_length:
                print("Warning: the dislocations are too close to select a new boundary (no more operation)!")
            else:
            
                """
                6:
                    Determine the middle of the biggest chunk,
                    which will tell us how to move the cells below and above 
                    the middle of the biggest chunk
                """
                mid_index = int(len(bin_biggest_chunk) / 2)
                mid_index = bin_biggest_chunk[mid_index]
                move_accumul = binList[mid_index] - z_low
                
                
                """
                7:
                    Move the cell below the middle interface up,
                    and move the cell above the middle interface down.
                    
                    Write the new cell in a file.
                """
                
                """
                #Adjust the particle coordinates by using OVITO (but sometimes it's prohibited')
                for i in range(len(data.particles.positions)):
                    if data.particles.positions[i,-1] <= z_high - move_accumul:
                        data.particles_.positions_[i,-1] = data.particles.positions[i,-1] + move_accumul
                    else:
                        data.particles_.positions[i,-1] = data.particles_.positions[i,-1] - (z_length - move_accumul)
                """
                
                final_coord = []
                
                
                #Only the z coordinates are changed        
                for i in range(len(data.particles.positions)):
                    final_coord.append([])
                    if data.particles.positions[i,-1] <= z_low + move_accumul:
                        final_coord[-1].append( data.particles.positions[i,0] )
                        final_coord[-1].append( data.particles.positions[i,1] )
                        final_coord[-1].append( data.particles.positions[i,2] + (z_length - move_accumul) )
                    else:
                        final_coord[-1].append( data.particles.positions[i,0] )
                        final_coord[-1].append( data.particles.positions[i,1] )
                        final_coord[-1].append( data.particles_.positions[i,-1] - move_accumul )
                                    
                  
                #Export a file by using OVITO (note the fist argument is "data" not the pipeline)
                #export_file(data, filein, "xyz", columns = ["Particle Type","Position.X", "Position.Y", "Position.Z"])
                
                totNum = data.particles.count
                
                lattice_text = "\"Lattice=\"%.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f\" Origin=\"%.6f %.6f %.6f\"" %\
                    (data.cell[0,0], data.cell[0,1], data.cell[0,2],\
                     data.cell[1,0], data.cell[1,1], data.cell[1,2],\
                     data.cell[2,0], data.cell[2,1], data.cell[2,2],\
                     data.cell[0,3], data.cell[1,3], data.cell[2,3])
                        
                # Access the property with the name 'Particle Type':
                prop = data.particles.particle_types
                
                with open(fileout, 'w') as out:
                    out.write("%d\n" % totNum)
                    out.write("%s\n" % lattice_text)
                    for i in range(len(final_coord)):
                        propid = prop[i]
                        out.write("%s %.5f %.5f %.5f %d\n" % \
                                  (prop.type_by_id(propid).name, final_coord[i][0], final_coord[i][1], final_coord[i][2], propid))
                        
                print("The write of a new cell is finished.")
                
def uncut_Ga2O3(filein, fileout, binWidth, limit_length, debug=False):
    """
    
    Similar to uncut, 
    but we need to remove Ga atoms,
    and only look at Oxygen sublattices which has FCC structures

    Parameters
    ----------
    filein : string
        The path to the input file.
    fileout : string
        The path to the output file
    binWidth : 
        The width of the bin that discritize the cell
    limit_length : 
        The minimum length of the region that is not covered by any dislocations.

    Returns
    -------
    None.
    (adjust the cell boundary)
    """
    
    # Enables activity logging for OVITO’s internal operations
    #ovito.enable_logging()
    
    pipeline = import_file(filein)

    data = pipeline.compute()
    
    """
    1:
      Discretize the cell along the z direction by a list of bin.
      
      We need to provide the width of the bin (may be modified by code),
      and a limit length from the function arguments.
      (If the uncovered region is smaller than this limit_length, we think
       the dislocations are too close and the cell boundaries will not be
       adjusted.)
    """
    
    #Get the lower limit and upper limit of the cell along the z direction
    z_low = data.cell[-1,-1]
    z_high = data.cell[-1,-1] + data.cell[-1,-2]
    z_length = data.cell[-1,-2]
    
    if debug == True:
        print(data.cell[:,:])
        
        print(data.particles.positions[:,:])
    
    binNumber = int(z_length / binWidth)
    
    #Recalculate the bin width according to the bin number
    binWidth = z_length / binNumber
    
    binList = np.arange(z_low, z_high, binWidth)
    
    """
    2:
        Remove Ga atoms
        
        Extract dislocations by using OVITO;
        Find dislocation segments.
        
        If no dislocation is cut by the z plane, flag=0, no operation;
        
        Otherwise, flag = 1, start the operation:
            split the dislocation cut by the z plane;
            move the newly generated dislocations according to the rule
            of periodic bondary condition.
            
        Note:
            if a dislocation is cut by the boundary, some dislocation segments
            calculated by OVITO will be out of the boundary limit.
            (The whole program is based on this property of OVITO)
        
    """
    
    #Select Ga atoms, whose types are 1 and 2
    pipeline.modifiers.append(
        SelectTypeModifier(
    operate_on = "particles",
    property = "Particle Type",
    types = { 1,2 }
        )
    )
    
    #Delete selected Ga atoms
    pipeline.modifiers.append( DeleteSelectedModifier() )
    
    # Extract dislocation lines from a crystal with FCC structure:
    modifier = DislocationAnalysisModifier()
    modifier.input_crystal_structure = DislocationAnalysisModifier.Lattice.FCC
    pipeline.modifiers.append(modifier)
    
    data = pipeline.compute()
    
    # Print list of dislocation lines:
    # Get the list of the z position of dislocation segments
    # (We only need the minimum and the maximum of the z position)
    flag = 0
    disloc_segment_z=[]
    #print("Found %i dislocation segments" % len(data.dislocations.segments))
    for segment in data.dislocations.segments:
        #print("Segment %i: length=%f, Burgers vector=%s" % (segment.id, segment.length, segment.true_burgers_vector))
        #print(segment.points)
        disz_mini = np.min(segment.points[:,-1])
        disz_max = np.max(segment.points[:,-1])
        
        if disz_mini >= z_low and disz_max <= z_high :
            disloc_segment_z.append([disz_mini, disz_max])
            
        elif disz_mini < z_low and disz_max > z_high :
            disloc_segment_z.append([z_low, z_high])
        
        #Split a dislocation to two, and apply the peridoic boundary condition    
        elif disz_mini < z_low and disz_max <= z_high :
            flag = 1
            disloc_segment_z.append([z_low, disz_max])
            disloc_segment_z.append([z_high - (z_low - disz_mini), z_high])
            
        #Split a dislocation to two, and apply the peridoic boundary condition    
        elif disz_mini >= z_low and disz_max > z_high :
            flag = 1
            disloc_segment_z.append([disz_mini, z_high])
            disloc_segment_z.append([z_low, z_low + (disz_max - z_high)])    
            
        else:
            print("Warning: special condition for disz_mini = %.5f and disz_max = %.5f." % (disz_mini, disz_max))
            print("In file: %s." % filein)
            continue
            
    if flag == 0:
        print("%s: no dislocation is cut by the z planes." % filein)
        #export xyz file
        pipeline = import_file(filein)
        data = pipeline.compute()
        totNum = data.particles.count
        lattice_text = "\"Lattice=\"%.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f\" Origin=\"%.6f %.6f %.6f\"" %\
                    (data.cell[0,0], data.cell[0,1], data.cell[0,2],\
                     data.cell[1,0], data.cell[1,1], data.cell[1,2],\
                     data.cell[2,0], data.cell[2,1], data.cell[2,2],\
                     data.cell[0,3], data.cell[1,3], data.cell[2,3])        
         # Access the property with the name 'Particle Type':
        prop = data.particles.particle_types
                
        with open(fileout, 'w') as out:
            out.write("%d\n" % totNum)
            out.write("%s\n" % lattice_text)
            for i in range(len(data.particles.positions)):
                propid = prop[i]
                out.write("%.5f %.5f %.5f %d\n" % \
                          (data.particles.positions[i][0], data.particles.positions[i][1], data.particles.positions[i][2], propid))  
        print("%s: transfer to xyz file." % fileout)     
                        
                        
               
    else:
        print("%s: Dislocation-cut is spotted." % filein)
        
        """
        3:
            Find the bins that are not covered by dislocation segments along the z direction,
            by looping over the list of dislocation segments (only containing the extreme segments).
            
            The index of the corresponding bins will be stored.
        """
        
        commen_bin_index = []
        for i in range(len(disloc_segment_z)):
            
            commen_bin_index.append([])
            
            for j in range(len(binList)):
                if not (binList[j] >= disloc_segment_z[i][0] and binList[j] <= disloc_segment_z[i][1]):
                    commen_bin_index[-1].append(j)
                    
        for i in range(len(commen_bin_index)):
            commen_bin_index[i] = np.array(commen_bin_index[i])
            
                
        """
        4:
            Compare the list of the bins that are not covered,
            in order to find the commen values.
            
            If no commen value is find, stop the operation;
            Otherwise, continue.
        """
        bin_noCover = commen_bin_index[0]
        
        for i in range(1, len(commen_bin_index)):
            bin_noCover = np.intersect1d(bin_noCover, commen_bin_index[i], assume_unique=True)
            
        
        if len(bin_noCover) < 1:
            print("Warning: all the space is covered by the dislocations! (no more operation)")
        else:
            
            """
            5:
                If the bins are connected with each other,
                they are organized to a chunk.
                
                If the biggest chunk is too small, then no more operation.
            """
            
            bin_chunks = []
            bin_chunks.append([])
            
            #Create the bin chunks
            if len(bin_noCover) == 1:
                bin_chunks[-1].append(bin_noCover[0])
            else:
                for i in range(1, len(bin_noCover)):
                    bin_chunks[-1].append(bin_noCover[i-1])
                    if bin_noCover[i] != (bin_noCover[i-1] + 1):
                        bin_chunks.append([])
                
                bin_chunks[-1].append(bin_noCover[-1])
                            
            
            #Find the biggest chunk
            bin_len_max = 0
            loc_max = 0
            for i in range(len(bin_chunks)):
                if len(bin_chunks[i]) > bin_len_max:
                    bin_len_max = len(bin_chunks[i])
                    loc_max = i
                    
            bin_biggest_chunk = bin_chunks[loc_max]
            
            if bin_len_max * binWidth < limit_length:
                print("Warning: the dislocations are too close to select a new boundary (no more operation)!")
            else:
            
                """
                6:
                    Determine the middle of the biggest chunk,
                    which will tell us how to move the cells below and above 
                    the middle of the biggest chunk
                """
                mid_index = int(len(bin_biggest_chunk) / 2)
                mid_index = bin_biggest_chunk[mid_index]
                move_accumul = binList[mid_index] - z_low
                
                
                """
                7:
                    Load again the cell
                    
                    Move the cell below the middle interface up,
                    and move the cell above the middle interface down.
                    
                    Write the new cell in a file.
                """
                
                """
                #Adjust the particle coordinates by using OVITO (but sometimes it's prohibited')
                for i in range(len(data.particles.positions)):
                    if data.particles.positions[i,-1] <= z_high - move_accumul:
                        data.particles_.positions_[i,-1] = data.particles.positions[i,-1] + move_accumul
                    else:
                        data.particles_.positions[i,-1] = data.particles_.positions[i,-1] - (z_length - move_accumul)
                """
                
                #Load again the cell, otherwise the Ga atoms are deleted.
                pipeline = import_file(filein)
                data = pipeline.compute()
                
                final_coord = []
                
                
                #Only the z coordinates are changed        
                for i in range(len(data.particles.positions)):
                    final_coord.append([])
                    if data.particles.positions[i,-1] <= z_low + move_accumul:
                        final_coord[-1].append( data.particles.positions[i,0] )
                        final_coord[-1].append( data.particles.positions[i,1] )
                        final_coord[-1].append( data.particles.positions[i,2] + (z_length - move_accumul) )
                    else:
                        final_coord[-1].append( data.particles.positions[i,0] )
                        final_coord[-1].append( data.particles.positions[i,1] )
                        final_coord[-1].append( data.particles_.positions[i,-1] - move_accumul )
                                    
                  
                #Export a file by using OVITO (note the fist argument is "data" not the pipeline)
                #export_file(data, filein, "xyz", columns = ["Particle Type","Position.X", "Position.Y", "Position.Z"])
                
                totNum = data.particles.count
                
                lattice_text = "\"Lattice=\"%.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f\" Origin=\"%.6f %.6f %.6f\"" %\
                    (data.cell[0,0], data.cell[0,1], data.cell[0,2],\
                     data.cell[1,0], data.cell[1,1], data.cell[1,2],\
                     data.cell[2,0], data.cell[2,1], data.cell[2,2],\
                     data.cell[0,3], data.cell[1,3], data.cell[2,3])
                        
                # Access the property with the name 'Particle Type':
                prop = data.particles.particle_types
                
                with open(fileout, 'w') as out:
                    out.write("%d\n" % totNum)
                    out.write("%s\n" % lattice_text)
                    for i in range(len(final_coord)):
                        propid = prop[i]
                        out.write("%.5f %.5f %.5f %d\n" % \
                                  (final_coord[i][0], final_coord[i][1], final_coord[i][2], propid))
                        
                print("The write of a new cell is finished.")

def modify_yz(filein,fileout):
    """change z and y axes
    filein:original dataframe
    fileout:z and y axes changed dataframe"""

    pipeline = import_file(filein)

    data = pipeline.compute()
    box1 = data.cell; #original cell
    
    data.cell_[2,2] = box1[1,1] #Interchange y and z
    data.cell_[1,1] = box1[2,2]
 
    smallbool = (abs(data.cell_[:]) < 1e-5) #optionally, set very small values of cell matrix to zero
    data.cell_[smallbool] = 0
    export_file(data, 'tempfile', "lammps/data",export_type_names=True )
    atom_number=81920
    # Open the input file in read mode
    with open('tempfile', 'r') as input_file:
    # Read the first 10 lines
        lines = input_file.readlines()[:17]

    # Open the output file in write mode
    with open(fileout, 'w') as output_file:
    # Write the first 10 lines into the output file
        output_file.writelines(lines)
    frame=pd.read_csv(filein, skiprows=17,sep=' ',nrows=atom_number,usecols=[0,1,2,3,4],header=None,names=['index','type', 'x', 'y', 'z'])
    temp_column = frame['y'].copy()
    frame['y']= frame['z']
    frame['z']=temp_column
    # Write the DataFrame to the CSV file without header
    frame.to_csv(fileout, sep=' ', index=False, mode='a', header=False)
  
    print(filein,"change yz cell is finished.")



if __name__ == '__main__':
    path='/home/heruhe/Desktop/Ga2o3/combine_cells/010_25boxes/boxes'
    os.chdir(path) 
    pka=[]
    for i in range(0,2000,10):
        if os.path.isfile('./data.lastframe-{}'.format(i)):
            pka.append(i)
            loc_Ga2O3='./data.lastframe-{}'.format(i)
            modify_frame = './frame-{}'.format(i)
            locout_Ga2O3='./{}.xyz'.format(i)
            modify_yz(loc_Ga2O3,modify_frame)
            uncut_Ga2O3(modify_frame, locout_Ga2O3, 1, 3.1652)
    print(pka)
    