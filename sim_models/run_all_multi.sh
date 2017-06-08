#!/bin/bash
#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=1:mem=100gb
#PBS -q pqtouldrid

## select Python script to run
RUNFILE=all_multi.py

## Copy from starting directory and move to space for running job:
cp $PBS_O_WORKDIR/$RUNFILE $TMPDIR 

## load Anaconda Python
module load anaconda3/4.3.1

## run the python file chosen
python $TMPDIR/$RUNFILE 

## copy outputted graphs back to original directory job script is located in
cp $TMPDIR/*.png $PBS_O_WORKDIR/all_multi_results