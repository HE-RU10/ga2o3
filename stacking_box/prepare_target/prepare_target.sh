#!/bin/bash

#Prepare target for the merging process
#1. Rotate: align the [010] direction with the z axis
#2. Prevent the cut of dislocations by the z boundary

echo "Preparing the target..."

mkdir Rotate

#My ovito module is in this virtual envrionment
#Please delete it
#source ~/Licument/virtualPython/spyder-env/bin/activate

echo "Rotating..."
python3 rotate.py

echo "Checking dislocations..."
python3 no_cut_dislocation.py > cutLog.txt

#Close the virtual environment
#deactivate 
