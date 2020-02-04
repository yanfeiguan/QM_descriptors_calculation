#!/bin/bash
#SBATCH -J qm_descriptors
#SBATCH -p medium
#SBATCH --time=24:00:00
#SBATCH --ntasks=40
#SBATCH --mem-per-cpu=2G
#SBATCH -N 1

PATH=$CONDA_PYTHON_PATH:$PATH:$NBOPATH
PYTHONPATH=$CONDA_PACKAGE_PATH:$PYTHONPATH

g16root=$G16_ROOT
GAUSS_SCRDIR=$SCRATH_FOLDER
export g16root GAUSS_SCRDIR
. $g16root/g16/bsd/g16.profile

python main.py
