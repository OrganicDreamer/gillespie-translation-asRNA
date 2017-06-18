#!/bin/bash
#PBS -l walltime=00:30:00
#PBS -l select=1:ncpus=1:mem=1gb

## select Python script to run
RUNFILE=gen-sweep-params1.py

## create/overwrite directory to store results
mkdir -p $WORK/sweep-params1

## Copy from starting directory and move to space for running job and generating results:
cp $PBS_O_WORKDIR/$RUNFILE $TMPDIR 
cp $PBS_O_WORKDIR/queued-production/1mR_finite_results/*.npz $TMPDIR

## load Anaconda Python
module load anaconda3/4.3.1

## run the python file chosen
python $TMPDIR/$RUNFILE 

## copy outputted data back to work directory results folder
cd $TMPDIR
cp ./*.npz $WORK/sweep-params1