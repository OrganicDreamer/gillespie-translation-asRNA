#!/bin/bash
#PBS -l walltime=03:00:00,mem=1gb

## load Anaconda Python

module load anaconda3/4.3.1

##select bash scripts to add to queue

qsub setup_run_py_script.sh
##qsub setup_run_py1_script.sh
##qsub setup_run_py2_script.sh
##etc..