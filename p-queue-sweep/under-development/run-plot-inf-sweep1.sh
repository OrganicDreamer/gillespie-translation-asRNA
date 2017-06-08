#!/bin/bash
#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=1:mem=100gb
#PBS -q pqtouldrid

## select Python script to run
RUNFILE=plot-inf-sweep1.py

## make directory to store outputted graphs
rm $PBS_O_WORKDIR/1mR-inf-sweep1-plots
mkdir $PBS_O_WORKDIR/1mR-inf-sweep1-plots

## Copy and move to space for running job:
cp $PBS_O_WORKDIR/$RUNFILE $TMPDIR 
cp $WORK/1mR-inf-sweep1-results $TMPDIR

## load Anaconda Python
module load anaconda3/4.3.1

## run the python file chosen
python $TMPDIR/$RUNFILE 

## copy outputted graphs back to original directory job script is located in
cp $TMPDIR/*.png $PBS_O_WORKDIR/1mR-inf-sweep1-plots

