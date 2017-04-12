#!/bin/bash
#PBS -l walltime=72:00:00,mem=100gb

## select Python script to run, copy it and move to work directory:
cp $HOME/1mR_finite.py $WORK ## select script to run here
cd $WORK

## load Anaconda Python and run script using Python
module load anaconda3/4.3.1
python 1mR_finite.py ## select script to run here

## copy outputted graphs back to home directory and erase from work directory
cp *.png $HOME
rm *.png