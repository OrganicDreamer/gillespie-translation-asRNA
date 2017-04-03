#!/bin/bash
#PBS -l walltime=03:00:00,mem=1gb
cp $HOME/1mR_finite.py $WORK 
cd $WORK
module load anaconda3/4.3.1
python 1mR_finite.py
cp *.png $HOME
rm *.png

