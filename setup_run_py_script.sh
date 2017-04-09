#!/bin/bash
#PBS -l walltime=03:00:00,mem=1gb

## select Python script to work on, copy it and move to work directory:
cp $HOME/1mR_finite.py $WORK ## select script to run here
cd $WORK

## run script using python (assuming correct module is already loaded):
python 1mR_finite.py ## select script to run here

## copy outputted graphs back to home directory and erase from work directory
cp *.png $HOME
rm *.png