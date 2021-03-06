#!/bin/bash
#PBS -l walltime=3:00:00
#PBS -l select=1:ncpus=1:mem=1gb
#PBS -q pqtouldrid

## select Python script to run
RUNFILE=1mR_finite.py

## Copy from starting directory and move to space for running job:
cp $PBS_O_WORKDIR/$RUNFILE $TMPDIR 

## load Anaconda Python
module load anaconda3/4.3.1

## run the python file chosen
python $TMPDIR/$RUNFILE 

## copy outputted graphs back to original directory job script is located in
cp $TMPDIR/*.png $PBS_O_WORKDIR/1mR_finite_results
cp $TMPDIR/*.npz $PBS_O_WORKDIR/1mR_finite_results