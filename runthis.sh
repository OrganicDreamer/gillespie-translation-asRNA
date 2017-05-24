#!/bin/bash
#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=1:mem=100gb

## select Python script to run, copy it and move to work directory:
cp $HOME/all_multi.py $WORK ## select script to run here
cd $WORK

## load Anaconda Python and run script using Python
module load anaconda3/4.3.1
python all_multi.py ## select script to run here

## copy outputted graphs back to home directory and erase from work directory
cp *.png $HOME
rm *.png