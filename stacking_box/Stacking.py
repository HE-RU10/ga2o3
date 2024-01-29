# -*- coding: utf-8 -*-
import os, tarfile
import numpy as np
import subprocess
from ovito.io import *
from ovito.modifiers import *
from ovito.data import DislocationNetwork
import ovito
import random
ovito.enable_logging()


x1 = 1
x2 = 71
k = 2


for i in range(x1,x2):
    try:
        os.remove(f'AT-{k}-{i}.PbcShrink.xyz')
    except:
        pass
    try:
        os.remove(f'AT-{k}-{i}.xyCorrected.xyz')
    except:
        pass
    try:
        os.remove(f'AT-{k}-{i}.DXACorrected.xyz')
    except:
        pass
    try:
        os.remove(f'AT-{k}-{i}.xyz')
    except:
        pass
try:
    os.remove(f'AT-{k}-stacking.xyz')
except:
    pass
try:
    os.remove('stacking.xyz')
except:
    pass



def untar(fname, dirs):

    try:
        t = tarfile.open(fname)
        t.extractall(path = dirs)
        return True
    except Exception as e:
        # print(e)
        pass
        return False

def PBC(x, y, z):
    Pos_x = np.nonzero(x > xborder)
    Neg_x = np.nonzero(x < -xborder)
    Pos_y = np.nonzero(y > yborder)
    Neg_y = np.nonzero(y < -yborder)       
    Pos_z = np.nonzero(z > zborder)
    Neg_z = np.nonzero(z < -zborder)
    x[Pos_x] = x[Pos_x] - 2.0 * xborder
    x[Neg_x] = x[Neg_x] + 2.0 * xborder
    y[Pos_y] = y[Pos_y] - 2.0 * yborder
    y[Neg_y] = y[Neg_y] + 2.0 * yborder
    z[Pos_z] = z[Pos_z] - 2.0 * zborder
    z[Neg_z] = z[Neg_z] + 2.0 * zborder
    return(x, y, z)

def mini_x_y(prog, ctype, file, eleCol, element, xDataCol, sizex, sizey, sizez):
    """
    Get the minimum x and y latice positions of a cell by using the FCC analyzer.
    Arguments: the arguments of fccAnalyzer
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
            print(f'Minimum lattice position: {min_x},{min_y},{min_z}')
            
        if "Reallocating" in lines[i]:
            continue
        elif "Allocating" in lines[i]:
            continue
        else:
            continue
    
    return(min_x, min_y, min_z)


def complement(wide, intervals):
    start, end = wide
    front,back = start, start
    result = []
    intervals = [sorted(x) for x in intervals]
    # print(intervals)
    # print(front,back)
    intervals.sort(key = lambda x: x[0])
    for (start1, end1) in intervals:
        # print(start1, end1)
        if start1 > end: break
        front = start1
        if front > back: result.append([back,front-0.1])
        if end1+0.1 > back: back=end1+0.1
    if back <= end: result.append([back,end])
    return(result)


prog = './minixy'
TotalNum = 0
label = 0
LatticeConstant = 3.1658
BoxxRef = 24*LatticeConstant*np.sqrt(6)
BoxyRef = 45*LatticeConstant*np.sqrt(2)
BoxzRef = 76*LatticeConstant*np.sqrt(3)/2

#-----------------------------------Scale+translation-----------------------------------

for i in range(x1, x2):
    try:
        untar(f'AT-{k}-{i}.xyz.tar.gz','./')
        x,y,z = np.loadtxt(f'AT-{k}-{i}.xyz',skiprows=2,usecols=[1,2,3],unpack=True)
        with open(f'AT-{k}-{i}.xyz', 'r+')as f1, open(f'AT-{k}-{i}.PbcShrink.xyz','a+') as f2, open(f'AT-{k}-{i}.xyCorrected.xyz','a+') as f3:
            lines = f1.readlines()
            print(lines[0], end = '', file = f2)
            print(lines[0], end = '', file = f3)
            TotalNum += int(lines[0])
            boxx = float(lines[1].split()[6])
            boxy = float(lines[1].split()[7])
            boxz = float(lines[1].split()[8])

            xborder = boxx/2.0
            yborder = boxy/2.0
            zborder = boxz/2.0

            xshift = random.uniform(1, 23)
            yshift = random.uniform(1, 44)
            x += xshift*LatticeConstant*np.sqrt(6)
            y += yshift*LatticeConstant*np.sqrt(2)

            # Random Mirroring
            ChangeOrNot_Y = random.choice([True, False])
            if ChangeOrNot_Y == True:
                y = boxy - y


            print(f'PBC modifying the frame {i}...')
            PBC(x, y, z)

            W = np.full(len(x),'W')
            # One = np.full(len(x),1)

            phix = BoxxRef/boxx
            phiy = BoxyRef/boxy
            phiz = BoxzRef/boxz

            print(f'frame {i},phix={phix},phiy={phiy},phiz={phiz}')

            x = x * phix + 0.3 #Eliminate the fluctuation of atomic layer after irradiation at the boundary, thus affecting subsequent link correction.
            y = y * phiy + 0.3
            z = z * phiz + 0.3
            
            # print(f'boxsize {BoxxRef} {BoxyRef} {BoxzRef}', file = f2)
            # print(f'boxsize {BoxxRef} {BoxyRef} {BoxzRef}', file = f3)
            print(f'Lattice=" {BoxxRef} 0.00000 0.00000 0.00000 {BoxyRef} 0.00000 0.00000 0.00000 {BoxzRef} " Origin=" {-0.5*BoxxRef} {-0.5*BoxyRef} {-0.5*BoxzRef} "', file = f2)
            print(f'Lattice=" {BoxxRef} 0.00000 0.00000 0.00000 {BoxyRef} 0.00000 0.00000 0.00000 {BoxzRef} " Origin=" {-0.5*BoxxRef} {-0.5*BoxyRef} {-0.5*BoxzRef} "', file = f3)

            datlstxyz = np.array([W,x,y,z]).T
            np.savetxt(f2, datlstxyz, fmt='%s',delimiter=' ', newline='\n')
            os.remove(f'AT-{k}-{i}.xyz')

        if i == x1:
            min_xRef, min_yRef, min_zCurrent = mini_x_y(prog, '-bcc', f'AT-{k}-{i}.PbcShrink.xyz', 1, 'W', 2, BoxxRef, BoxyRef, BoxzRef)
            min_xCurrent, min_yCurrent= min_xRef, min_yRef

        else:
            min_xCurrent, min_yCurrent, min_zCurrent = mini_x_y(prog, '-bcc', f'AT-{k}-{i}.PbcShrink.xyz', 1, 'W', 2, BoxxRef, BoxyRef, BoxzRef)

        deltaX = min_xCurrent - min_xRef
        deltaY = min_yCurrent - min_yRef
        print(f'frame{i},deltaX={deltaX},deltaY={deltaY}')

        x,y,z = np.loadtxt(f'AT-{k}-{i}.PbcShrink.xyz',delimiter=' ',skiprows=2,usecols=[1,2,3],unpack=True)
        os.remove(f'AT-{k}-{i}.PbcShrink.xyz')

        x -= deltaX
        y -= deltaY
        xborder = BoxxRef/2.0
        yborder = BoxyRef/2.0
        zborder = BoxzRef/2.0

        print(f'PBC modifying the frame {i} again...')
        PBC(x, y, z)
        datlstxyz = np.array([W,x,y,z]).T
        with open(f'AT-{k}-{i}.xyCorrected.xyz','a+') as f1:
            np.savetxt(f1, datlstxyz, fmt='%s',delimiter=' ', newline='\n')
        print(f'-----------------frame {i} XY coreection, done-----------------')
        if i > x1:
            label += 1
        
    except:
        pass



#-----------------------------------DXA judgement-----------------------------------
print('*'*100)
print(f'DXA judgement started...')
modifier = DislocationAnalysisModifier()
modifier.input_crystal_structure = DislocationAnalysisModifier.Lattice.BCC
zborder = BoxzRef/2

for i in range(x1, x2):
    # try:
            pipeline = import_file(f'AT-{k}-{i}.xyCorrected.xyz')
            pipeline.modifiers.append(modifier)
            data = pipeline.compute()
            x,y,z = np.loadtxt(f'AT-{k}-{i}.xyCorrected.xyz',delimiter=' ',skiprows=2,usecols=[1,2,3],unpack=True)
            os.remove(f'AT-{k}-{i}.xyCorrected.xyz')

            RangePerDis = []


            print("Found %i dislocation segments" % len(data.dislocations.segments))
            for segment in data.dislocations.segments:
                # print("Segment %i: length=%f, Burgers vector=%s" % (segment.id, segment.length, segment.true_burgers_vector))
               
               SegZMin = (segment.points[:,-1]).min()
               SegZMax = (segment.points[:,-1]).max() 
               # SegZMin and SegZMax respresent the zlo and zhi positions of each dislocation
               if SegZMin < -zborder:
                  deltaz = abs(SegZMin) - zborder
                  sep1 = [-zborder, SegZMax]
                  sep2 = [zborder-deltaz, zborder]
                  RangePerDis.append(sep1)
                  RangePerDis.append(sep2)
               elif SegZMax > zborder:
                  deltaz = SegZMax - zborder
                  sep1 = [-zborder, -zborder+deltaz]
                  sep2 = [SegZMin, zborder]
                  RangePerDis.append(sep1)
                  RangePerDis.append(sep2)
               else:
                  perdis = [SegZMin, SegZMax]
                  RangePerDis.append(perdis)


            ZRange = [-zborder, zborder]
            ComplementRange = complement(ZRange, RangePerDis)
            RandomPickIndex = round(len(ComplementRange)-1)
            if RandomPickIndex < 0:
                print(f'For frame {i}, the cell is no need to be shifted.')
            else:
                widthlst = []
                for p in range(len(ComplementRange)):
                    Width = ComplementRange[p][1]-ComplementRange[p][0]
                    widthlst.append(Width)
                WidestIndex = np.argmax(widthlst)

                PickPositionZlo = ComplementRange[WidestIndex][0] + np.max(widthlst)/5.0
                PickPositionZhi = ComplementRange[WidestIndex][1] - np.max(widthlst)/5.0

                shiftPosition = random.uniform(PickPositionZlo, PickPositionZhi)
                shiftValue = zborder - shiftPosition
                z += shiftValue 
                print(f'Pick shift position from range({PickPositionZlo}, {PickPositionZhi}).')
                print(f'For frame {i}, boundary thickness={np.max(widthlst)/5.0}.')
                print(f'The cell has been shifted towards +z with a value of {shiftValue}.')
                widthlst = []

            PBC(x, y, z)

            W = np.full(len(x),'W')
            datlstxyz = np.array([W,x,y,z]).T
            with open(f'AT-{k}-{i}.DXACorrected.xyz','a+') as f1:
                print(len(x), file = f1)
                print(f'Lattice=" {BoxxRef} 0.00000 0.00000 0.00000 {BoxyRef} 0.00000 0.00000 0.00000 {BoxzRef} " Origin=" {-0.5*BoxxRef} {-0.5*BoxyRef} {-0.5*BoxzRef} "', file = f1)
                np.savetxt(f1, datlstxyz, fmt='%s',delimiter=' ', newline='\n')

            x,y,z = np.loadtxt(f'AT-{k}-{i}.DXACorrected.xyz',delimiter=' ',skiprows=2,usecols=[1,2,3],unpack=True)
            

            if i == x1:
                min_xRef, min_yRef, min_zCurrent = mini_x_y(prog, '-bcc', f'AT-{k}-{i}.DXACorrected.xyz', 1, 'W', 2, BoxxRef, BoxyRef, BoxzRef)
                min_xCurrent, min_yCurrent= min_xRef, min_yRef
            else:
                min_xCurrent, min_yCurrent, min_zCurrent = mini_x_y(prog, '-bcc', f'AT-{k}-{i}.DXACorrected.xyz', 1, 'W', 2, BoxxRef, BoxyRef, BoxzRef)
            os.remove(f'AT-{k}-{i}.DXACorrected.xyz')
            deltaZ = min_zCurrent + zborder
            print(f'frame{i},deltaZ={deltaZ}')
            z -= deltaZ + 0.4

            PBC(x, y, z)

            datlstxyz = np.array([W,x,y,z]).T
            with open(f'AT-{k}-{i}.DXACorrected.xyz','a+') as f1:
                print(len(x), file = f1)
                print(f'Lattice=" {BoxxRef} 0.00000 0.00000 0.00000 {BoxyRef} 0.00000 0.00000 0.00000 {BoxzRef} " Origin=" {-0.5*BoxxRef} {-0.5*BoxyRef} {-0.5*BoxzRef} "', file = f1)
                np.savetxt(f1, datlstxyz, fmt='%s',delimiter=' ', newline='\n')

            # Correct the lattice structure link at Z border

            min_xCurrent, min_yCurrent, min_zCurrent = mini_x_y(prog, '-bcc', f'AT-{k}-{i}.DXACorrected.xyz', 1, 'W', 2, BoxxRef, BoxyRef, BoxzRef)
            x,y,z = np.loadtxt(f'AT-{k}-{i}.DXACorrected.xyz',delimiter=' ',skiprows=2,usecols=[1,2,3],unpack=True)
            os.remove(f'AT-{k}-{i}.DXACorrected.xyz')

            Min_Z_layer = np.nonzero(z<(min_zCurrent+0.5))
            Min_Z_layer_index = np.argwhere(z<(min_zCurrent+0.5))
            X_on_minZlst = np.zeros(0)

            for index in range(len(Min_Z_layer_index)):
                X_on_minZ = x[Min_Z_layer_index[index]]
                X_on_minZlst = np.append(X_on_minZlst, X_on_minZ)

            minX_on_minZ_index = np.nonzero(X_on_minZlst < X_on_minZlst.min()+3)
            ave_x_on_minz = np.mean(X_on_minZlst[minX_on_minZ_index])

            print(f'min_xCurrent={min_xCurrent},ave_x_on_minz={ave_x_on_minz}')

            if abs(abs(min_xCurrent)-abs(ave_x_on_minz)) > 1.8:
                print(f'frame {i} needs to be corrected 0.87*2 for the bottom link...')
                z -= 0.87*2 #indicate the min X on min Z layer is body center atom
            elif abs(abs(min_xCurrent)-abs(ave_x_on_minz)) > 0.6:
                print(f'frame {i} needs to be corrected 0.87*2 for the bottom link...')
                z -= 0.87 #indicate the min X on min Z layer is body center atom
            PBC(x, y, z)
            datlstxyz = np.array([W,x,y,z]).T

            with open(f'AT-{k}-{i}.DXACorrected.xyz','a+') as f1:
                print(len(x), file = f1)
                print(f'Lattice=" {BoxxRef} 0.00000 0.00000 0.00000 {BoxyRef} 0.00000 0.00000 0.00000 {BoxzRef} " Origin=" {-0.5*BoxxRef} {-0.5*BoxyRef} {-0.5*BoxzRef} "', file = f1)
                np.savetxt(f1, datlstxyz, fmt='%s',delimiter=' ', newline='\n')




            if i == x1:
                min_xRef, min_yRef, min_zCurrent = mini_x_y(prog, '-bcc', f'AT-{k}-{i}.DXACorrected.xyz', 1, 'W', 2, BoxxRef, BoxyRef, BoxzRef)
                min_yCurrent = min_yRef
                min_xCurrent = min_xRef

            else:
                min_xCurrent, min_yCurrent, min_zCurrent = mini_x_y(prog, '-bcc', f'AT-{k}-{i}.DXACorrected.xyz', 1, 'W', 2, BoxxRef, BoxyRef, BoxzRef)

            deltaX = min_xCurrent - min_xRef
            deltaY = min_yCurrent - min_yRef

            print(f'frame{i},deltaX={deltaX},deltaY={deltaY}')

            x,y,z = np.loadtxt(f'AT-{k}-{i}.DXACorrected.xyz',delimiter=' ',skiprows=2,usecols=[1,2,3],unpack=True)
            os.remove(f'AT-{k}-{i}.DXACorrected.xyz')

            x -= deltaX
            y -= deltaY

            xborder = BoxxRef/2.0
            yborder = BoxyRef/2.0
            zborder = BoxzRef/2.0

            print(f'PBC modifying the frame {i} again...')
            PBC(x, y, z)
            datlstxyz = np.array([W,x,y,z]).T
            with open(f'AT-{k}-{i}.DXACorrected.xyz','a+') as f1:
                print(len(x), file = f1)
                print(f'Lattice=" {BoxxRef} 0.00000 0.00000 0.00000 {BoxyRef} 0.00000 0.00000 0.00000 {BoxzRef} " Origin=" {-0.5*BoxxRef} {-0.5*BoxyRef} {-0.5*BoxzRef} "', file = f1)
                np.savetxt(f1, datlstxyz, fmt='%s',delimiter=' ', newline='\n')



            print(f'-----------------frame {i} DXA+link coreection, done-----------------')

    # except:
        # pass 



#-----------------------------------Stacking-----------------------------------
print('*'*100)
boxz = 0
gapz = 0

print('Writing the final stacking cell >.< ')
for i in range(x1, x2):
    try:
        x,y,z = np.loadtxt(f'AT-{k}-{i}.DXACorrected.xyz',delimiter=' ',skiprows=2,usecols=[1,2,3],unpack=True)
        with open(f'AT-{k}-{i}.DXACorrected.xyz','r+') as f1, open(f'AT-{k}-stacking.xyz', 'a+')as f2:
            lines = f1.readlines()
            if i == x1:
                half = float(lines[1].split()[9])/2
                boxz = half
                z -= half
            if i > x1:
                boxz += float(lines[1].split()[9]) + gapz
                z -= boxz

            datlstdata = np.array([x,y,z]).T
            W = np.full(len(x),'W').reshape((len(x), 1))
            One = np.full(len(x),1).reshape((len(x), 1))

            datlstAll = np.concatenate((W, datlstdata, One), axis=1)
            np.savetxt(f2, datlstAll, fmt='%s',delimiter=' ', newline='\n')
    except:
        pass


with open(f'AT-{k}-stacking.xyz','r+') as f1:
    content = f1.read()
    f1.seek(0, 0)
    f1.write(f'Lattice=" {BoxxRef} 0.00000 0.00000 0.00000 {BoxyRef} 0.00000 0.00000 0.00000 {boxz+half} " Origin=" {-0.5*BoxxRef} {-0.5*BoxyRef} {-boxz-half} "\n'+content)

with open(f'AT-{k}-stacking.xyz','r+') as f1:
    content = f1.read()
    f1.seek(0, 0) 
    f1.write(f'{TotalNum}\n'+content)



for k in range(1, 4):
    for i in range(x1, x2):
        try:
            os.remove(f'AT-{k}-{i}.DXACorrected.xyz')
        except:
            pass






