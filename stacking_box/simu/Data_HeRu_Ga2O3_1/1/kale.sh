#!/bin/bash
#SBATCH --job-name=Ga2O3_1
#SBATCH -o screen1.txt              ###Output file
#SBATCH -c 1                        ###Core count
#SBATCH -n 50                       ###Number of tasks
#SBATCH -N 1                       ###Number of nodes
#SBATCH -t 2-23:59:59                 ###Requested time
#SBATCH --mem=50G                    ###Requested memory
#SBATCH --mail-type=END
#SBATCH --mail-user=xin.jin@helsinki.fi

module purge                                ## Purge modules for a clean start
#module load Python/3.7.0-intel-2018b
module load Python/3.10.4-GCCcore-11.3.0
#module load OpenMPI/4.0.1-GCC-8.3.0-2.32 #Can get segmentation fault
module load OpenMPI/4.0.3-GCC-9.3.0

source /home/jinxinjx/proj/pyVirtual/bin/activate

python3 merRelax.py
python3 quench_rbs.py
python3 rbs_simu1.py

deactivate

