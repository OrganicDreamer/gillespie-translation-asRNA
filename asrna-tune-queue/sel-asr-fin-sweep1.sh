#!/bin/bash
#PBS -l walltime=01:00:00
#PBS -l select=1:ncpus=1:mem=1gb
#PBS -J 0-24
#PBS -o myout
#PBS -e myerr
#PBS -j oe


## select Python script to run
RUNFILE=sel-asr-fin-sweep1.py

## create/overwrite directory to store results
mkdir -p $PBS_O_WORKDIR/sel-asr-fin-sweep1

## Copy and move to space for running job:
cp $PBS_O_WORKDIR/$RUNFILE $TMPDIR 
cp $WORK/asr-fin-sweep-results1/*.npz $TMPDIR

## load Anaconda Python
module load anaconda3/4.3.1

## run the python file chosen
python $TMPDIR/$RUNFILE 

## copy outputted graphs back to original directory job script is located in
cd $TMPDIR
cp ./*.png $PBS_O_WORKDIR/sel-asr-fin-sweep1
cp ./*asr-on-table-entry.npz $PBS_O_WORKDIR/sel-asr-fin-sweep1
cp ./fin-asrna.npz $PBS_O_WORKDIR/sel-asr-fin-sweep1
