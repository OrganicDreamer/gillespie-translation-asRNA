#!/bin/bash
#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=1:mem=24gb
#PBS -J 0-49

## select Python script to run
RUNFILE=1mr-fin-sweep.py

## create/overwrite directory to store results
mkdir -p $WORK/1mr-fin-sweep-results

## Copy from starting directory and move to space for running job and generating results:
cp $PBS_O_WORKDIR/$RUNFILE $TMPDIR 

## load Anaconda Python
module load anaconda3/4.3.1

## run the python file chosen
python $TMPDIR/$RUNFILE 

## copy outputted data back to work directory results folder
cd $TMPDIR
cp ./*.npz $WORK/1mr-fin-sweep-results