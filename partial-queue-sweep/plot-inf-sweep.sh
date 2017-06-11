#!/bin/bash
#PBS -l walltime=01:00:00
#PBS -l select=1:ncpus=1:mem=1gb

## select Python script to run
RUNFILE=plot-inf-sweep.py

## create/overwrite directory to store results
mkdir -p $PBS_O_WORKDIR/1mr-inf-sweep-plots

## Copy and move to space for running job:
cp $PBS_O_WORKDIR/$RUNFILE $TMPDIR 
cp $WORK/1mr-inf-sweep-results/*.npz $TMPDIR

## load Anaconda Python
module load anaconda3/4.3.1

## run the python file chosen
python $TMPDIR/$RUNFILE 

## copy outputted graphs back to original directory job script is located in
cd $TMPDIR
cp ./*.png $PBS_O_WORKDIR/1mr-inf-sweep-plots

