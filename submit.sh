#!/bin/bash
#$ -N test
#$ -l long
#$ -l h_rt=120:00:00
##$ -l harpertown
#$ -m ae
#$ -pe singlenode  24
#$ -cwd
#$ -o out.txt
#$ -e err.txt

PATH=$CONDA_PYTHON_PATH:$PATH:/home/yanfeig:$NBOPATH
PYTHONPATH=$CONDA_PACKAGE_PATH:$PYTHONPATH

g16root=$G16ROOT
GAUSS_SCRDIR=$SCRATCH_FOLDER_G16
export g16root GAUSS_SCRDIR
. $g16root/g16/bsd/g16.profile

python main.py --ismiles test.csv

