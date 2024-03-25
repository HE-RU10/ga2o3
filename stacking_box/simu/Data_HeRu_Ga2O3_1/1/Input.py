#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  5 13:01:31 2019

@author: jin

Write the RBSADEC input file

Path:../Modules/RBSout/Input.py

Last modified: 2019-11-13
"""

class input_dat():
    """
    RBSADEC input.dat
    """
    
    def __init__(self):
        
        self.flag = '-rbs'
        
        self.Ion = 'He'
        self.Energy = 1e6 #[eV]
        self.Nhistory = 1e3
        
        
        self.Theta_min = 0.0
        self.Theta_max = 0.1
        self.Fii_min = 0.0
        self.Fii_max = 360.0
        self.Rotation_theta = 0.0
        self.Rotation_fii = 0.0
        self.Dechannel_theta = 0.0
        self.Mode_rotate = 0
        self.Mode_ibd = 1
        self.Mode_tbm = 0
        
        
        self.Ntype = 1
        self.Density = 2.329 #[g/cm3]
        self.Element = ['Si']
        self.Displacement = [0.10] #[\AA]
        self.Mode_vibrate = 1
        self.DebyeDisp = 1
        self.DebyeTemp = 300 #[K]
        
        
        self.Mini_potential = 0 #[eV]
        self.Coll_offset = 0 #[\AA]
        self.Free_distance = 100 #[\AA]
        self.Elec_distance = 100 #[\AA]
        self.Threshold_energy = 50 #[eV]
        self.Mode_FFP = 0
        self.Mode_Impact = 1
        self.Elstop_mode = 0
        
        
        self.Positionmin_x = 0 #[\AA]
        self.Positionmin_y = 0 #[\AA]
        self.Positionmin_z = -3 #[\AA]
        self.Positionmax_x = 100 #[\AA]
        self.Positionmax_y = 100 #[\AA]
        self.Positionmax_z = -3 #[\AA]
        
        self.xlow = 0 #[\AA]
        self.ylow = 0 #[\AA]
        self.zlow = 0 #[\AA]
        self.Boxsize_x = 100 #[\AA]
        self.Boxsize_y = 100 #[\AA]
        self.Boxsize_z = 100 #[\AA]
        
        self.Depthsize = 1e3 #[\AA]
        self.DepthStop = 1e3 #[\AA]
        
        self.Boundary_x = 1
        self.Boundary_y = 1
        self.Boundary_z = 0
        
        self.Binsize_x = 5 #[\AA]
        self.Binsize_y = 5 #[\AA]
        self.Binsize_z = 5 #[\AA]
        
        self.DBinsize_x = 5 #[\AA]
        self.DBinsize_y = 5 #[\AA]
        self.DBinsize_z = 5 #[\AA] 
        
        
        self.Temperature = 285 #[K]
        self.SpreadAngle = 60
        self.Seed = 2579
        
        
        self.Detector_theta = 165
        self.Detector_fii = 0
        self.Detector_radius = 10 #[mm]
        self.Detector_distance = 5 #[cm]
        self.Nchannel = 600
        self.Detector_FWHM = 15 #[keV]
        self.Detector_slope = 2 #[keV/channel]
        self.Detector_intercept = 0 #[keV]
        
        
        self.table_dat = './Data/table.dat'
        self.scatter_in = './Data/Scatter.mat'
        self.time_in = './Data/Time.mat'
        self.elstop1_in = './in/elstop1.in'
        self.coords1_in = './in/coords1.in'
        
        self.nbox = 1        
        self.regionBox = []
        self.nsphere = 1
        self.regionSphere = []
        
        self.NRA_mode = 0
        self.NRA_nucleus = 'O16'
        self.NRA_emit = 'He4'
        self.NRA_residual = 'N14'
        self.NRA_Q = 1e6
        self.NRA_spreadAngle = 0.0
        self.NRA_crossSection = './in/diffCS.in'
        self.NRA_elstop = './in/elstop_nra.in'
        
    def write_dat(self, loc):
        """
        Write the input.dat for RBSADEC
        """
        
        with open(loc, 'w') as wt:
            
            wt.write("*****************************************************\n")
            wt.write("*                   RBSADEC INPUT                   *\n")
            wt.write("*****************************************************\n")
            wt.write("\n")
            wt.write("Ion             =   %s\n"%(self.Ion))
            wt.write("Energy          =   %.1f #[eV]\n"%(self.Energy))
            wt.write("Nhistory        =   %d\n"%(self.Nhistory))
            wt.write("Theta.min       =   %.2f\n"%(self.Theta_min))
            wt.write("Theta.max       =   %.2f #IBD characteristic angle\n"%(self.Theta_max))
            wt.write("Fii.min         =   %.2f\n"%(self.Fii_min))
            wt.write("Fii.max         =   %.2f\n"%(self.Fii_max))
            wt.write("Rotation.theta  =   %.2f\n"%(self.Rotation_theta))
            wt.write("Rotation.fii    =   %.2f\n"%(self.Rotation_fii))
            wt.write("Dechannel.theta =   %.2f\n"%(self.Dechannel_theta))
            wt.write("Mode.rotate     =   %d #0: go to Rotation.fii; 1: fii from 0 to 359; 2: based on 1, detector also rotates\n"%(self.Mode_rotate))
            wt.write("Mode.ibd        =   %d #0: Type A (original); 1: Type B; 2: Type C\n"%(self.Mode_ibd))
            wt.write("Mode.tbm        =   %d #0: no action; 1: distinguish the channeled and dechanneled ion\n"%(self.Mode_tbm))
            wt.write("*****************************************************\n")
            wt.write("\n")
            wt.write("Ntype            =  %d\n"%(self.Ntype))
            wt.write("Density          =  %.3f #[g/cm3]\n"%(self.Density))
            for i in range(len(self.Element)):
                wt.write("Element.%d        =  %s\n"%(i+1, self.Element[i]))
            
            wt.write('DebyeDisp        =  %d #Recommend 0 for MD cells\n'%(self.DebyeDisp))
            wt.write("Mode.vibrate     =  %d #0 & 2: calculated by Debye theory; 1: read from inputs\n"%(self.Mode_vibrate))
            if(1 == self.Mode_vibrate):
                for i in range(len(self.Displacement)):
                    wt.write("Displacement.%d   =  %.3f\n"%(i+1, self.Displacement[i]))
            elif (0 == self.Mode_vibrate or 2 == self.Mode_vibrate):
                wt.write("DebyeT           = %.1f #[K]\n"%(self.DebyeTemp))
                wt.write("Temperature      =  %.1f #[K]\n"%(self.Temperature))
            else:
                raise ValueError("The value of Mode.vibrate can only be 0, 1 or 2.")
                
            wt.write("\n")
            wt.write("Mini.potential   =  %.3f #[eV]\n"%(self.Mini_potential))
            wt.write("Coll.offset      =  %.3f\n"%(self.Coll_offset))
            wt.write("Free.distance    =  %.1f\n"%(self.Free_distance))
            wt.write("Threshold.energy =  %.1f\n"%(self.Threshold_energy))
            wt.write("Mode.FFP         =  %d #0: Poisson process (original); 1: fixed distance; 2: TRIM\n"%(self.Mode_FFP))
            wt.write("Mode.Impact      =  %d #0: Iradina; 1: TRIM; 2: Original\n"%(self.Mode_Impact))
            wt.write("Elstop.mode      =  %d\n"%(self.Elstop_mode))
            wt.write("*****************************************************\n")
            wt.write("\n")
            wt.write("Positionmin.x   =  %.5f\n"%(self.Positionmin_x))
            wt.write("Positionmin.y   =  %.5f\n"%(self.Positionmin_y))
            wt.write("Positionmin.z   =  %.5f\n"%(self.Positionmin_z))
            wt.write("\n")
            wt.write("Positionmax.x   =  %.5f\n"%(self.Positionmax_x))
            wt.write("Positionmax.y   =  %.5f\n"%(self.Positionmax_y))
            wt.write("Positionmax.z   =  %.5f\n"%(self.Positionmax_z))
            wt.write("\n")
            wt.write("Box.xlow        =  %.5f\n"%(self.xlow))
            wt.write("Box.ylow        =  %.5f\n"%(self.ylow))
            wt.write("Box.zlow        =  %.5f\n"%(self.zlow))
            wt.write("Boxsize.x       =  %.5f\n"%(self.Boxsize_x))
            wt.write("Boxsize.y       =  %.5f\n"%(self.Boxsize_y))
            wt.write("Boxsize.z       =  %.5f\n"%(self.Boxsize_z))
            wt.write("\n")
            wt.write("Depthsize       =  %.5f\n"%(self.Depthsize))
            wt.write("DepthStop       =  %.5f\n"%(self.DepthStop))
            wt.write("\n")
            wt.write("Boundary.x      =  %d\n"%(self.Boundary_x))
            wt.write("Boundary.y      =  %d\n"%(self.Boundary_y))
            wt.write("Boundary.z      =  %d\n"%(self.Boundary_z))
            wt.write("\n")
            wt.write("Binsize.x       =  %.5f\n"%(self.Binsize_x))
            wt.write("Binsize.y       =  %.5f\n"%(self.Binsize_y))
            wt.write("Binsize.z       =  %.5f\n"%(self.Binsize_z))
            wt.write("\n")
            if('-loc' in self.flag):
                wt.write("DBinsize.x      =  %.5f\n"%(self.DBinsize_x))
                wt.write("DBinsize.y      =  %.5f\n"%(self.DBinsize_y))
                wt.write("DBinsize.z      =  %.5f\n"%(self.DBinsize_z))
            wt.write("*****************************************************\n")
            wt.write("\n")
            wt.write("SpreadAngle     =  %.1f\n"%(self.SpreadAngle))
            wt.write("Seed            =  %d\n"%(self.Seed))
            wt.write("*****************************************************\n")
            wt.write("\n")
            wt.write("Detector.theta     =  %.2f\n"%(self.Detector_theta))
            wt.write("Detector.fii       =  %.2f\n"%(self.Detector_fii))
            wt.write("\n")
            wt.write("Detector.radius    =  %.1f #[mm]\n"%(self.Detector_radius))
            wt.write("Detector.distance  =  %.1f #[cm]\n"%(self.Detector_distance))
            wt.write("Detector.FWHM      =  %.1f #[keV]\n"%(self.Detector_FWHM))
            wt.write("\n")
            wt.write("Nchannel           =  %d\n"%(self.Nchannel))
            wt.write("Detector.slope     =  %.3f #[keV/channel]\n"%(self.Detector_slope))
            wt.write("Detector.intercept =  %.3f #[keV]\n"%(self.Detector_intercept))
            wt.write("*****************************************************\n")
            wt.write("\n")
            wt.write("table.dat    = %s\n"%(self.table_dat))
            wt.write("scatter.in   = %s\n"%(self.scatter_in))
            wt.write("time.in      = %s\n"%(self.time_in))
            wt.write("\n")
            wt.write("elstop1.in   = %s\n"%(self.elstop1_in))
            wt.write("coords1.in   = %s\n"%(self.coords1_in))
            wt.write("*****************************************************\n")
            wt.write("\n")
            if('-loc2' in self.flag):
                wt.write("Nbox   = %d\n"%(self.nbox))
                for j in range(len(self.regionBox)):
                    wt.write("Region.box.%d =" % (j+1) )
                    for i in range(len(self.regionBox[j])):
                        wt.write(" %.3f" % (self.regionBox[j][i]))
                    wt.write("\n")
                
                wt.write("Nsphere   = %d\n"%(self.nsphere))
                for j in range(len(self.regionSphere)):
                    wt.write("Region.sph.%d =" % (j+1) )
                    for i in range(len(self.regionSphere[j])):
                        wt.write(" %.3f" % (self.regionSphere[j][i]))
                    wt.write("\n")
                
                wt.write("*****************************************************\n")
                wt.write("\n")
                    
            
            if ('-nra' in self.flag) :
                wt.write("NRA.mode       =  %d #0: Resonance elastic scattering; 1: Nuc. react.\n"%(self.NRA_mode))
                wt.write("NRA.nucleus    =  %s\n"%(self.NRA_nucleus))
                wt.write("NRA.emit       =  %s\n"%(self.NRA_emit))
                wt.write("NRA.residual   =  %s\n"%(self.NRA_residual))
                wt.write("NRA.Q          =  %s\n"%(self.NRA_Q))
                wt.write("NRA.sa         =  %.3f\n"%(self.NRA_spreadAngle))
                wt.write("NRA.csection   =  %s\n"%(self.NRA_crossSection))
                wt.write("NRA.elstop     =  %s\n"%(self.NRA_elstop))
                wt.write("*****************************************************\n")
                wt.write("\n")

if __name__ == '__main__':
    
    """Test 2019-11-05"""
    TestA = input_dat()
    TestA.write_dat('input.dat')
