#!/bin/bash
#PBS -l walltime=03:00:00,mem=1gb

## select Python script to work on, copy it and move to work directory:

cp $HOME/1mR_finite.py $WORK ## modify script to run here
cd $WORK

## set up and load python and interpreter and run script:

module load anaconda3/4.3.1
python 1mR_finite.py ## modify script to run here

## copy outputted graphs back to home directory and erase from work directory
cp *.png $HOME
rm *.png

