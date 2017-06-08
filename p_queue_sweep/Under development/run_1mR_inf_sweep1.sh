#!/bin/bash
#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=1:mem=100gb
#PBS -q pqtouldrid
#PBS -J 0-49

## select Python script to run
RUNFILE=1mR_inf_sweep1.py

## Copy from starting directory and move to space for running job and storing results:
cp $PBS_O_WORKDIR/$RUNFILE $TMPDIR 

## load Anaconda Python
module load anaconda3/4.3.1

## run the python file chosen
python $TMPDIR/$RUNFILE 

## copy outputted data back to work directory results folder
cp $TMPDIR/*.npz $WORK/1mR_inf_sweep1_results