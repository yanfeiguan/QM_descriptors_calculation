#!/bin/bash
#$ -N Ezlopitant
#$ -l long
#$ -l h_rt=120:00:00
#$ -l h="node70|node71|node72|node78|node78|node80|node82|node83|node84|node88|node89"
##$ -l harpertown
#$ -m ae
#$ -pe singlenode  24
#$ -cwd
#$ -o out.txt
#$ -e err.txt

PATH=/home/yanfeig/anaconda3/envs/general/bin:$PATH:/home/yanfeig:/home/yanfeig/nbo6/bin
PYTHONPATH=/home/yanfeig/anaconda3/envs/general/lib/python3.5/site-packages:$PYTHONPATH

g16root=/opt
GAUSS_SCRDIR=/scratch/yanfeig
export g16root GAUSS_SCRDIR
. $g16root/g16/bsd/g16.profile

xtbroot=/home/yanfeig/software/xtb_exe/bin

python main.py

