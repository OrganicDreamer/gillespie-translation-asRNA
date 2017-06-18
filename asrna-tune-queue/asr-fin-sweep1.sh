#!/bin/bash
#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=1:mem=12gb
#PBS -J 0-624
#PBS -o myout
#PBS -e myerr
#PBS -j oe

## select Python script to run
RUNFILE=asr-fin-sweep1.py

## create/overwrite directory to store results
mkdir -p $WORK/asr-fin-sweep-results1

## Copy from starting directory and move to space for running job and generating results:
cp $PBS_O_WORKDIR/$RUNFILE $TMPDIR 
cp $WORK/sweep-params1/*.npz $TMPDIR

## load Anaconda Python
module load anaconda3/4.3.1

## run the python file chosen
python $TMPDIR/$RUNFILE 

## copy outputted data back to work directory results folder
cd $TMPDIR
cp ./*fin-sweep.npz $WORK/asr-fin-sweep-results1
cp ./fin-asrna.npz $WORK/asr-fin-sweep-results1