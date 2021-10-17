#!/bin/bash
#SBATCH -A lilyw7
#SBATCH --partition=free
#SBATCH -t 24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=48
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4gb
#SBATCH --requeue

## USAGE

###sbatch -J <JOBNAME> -o <stdoutfile.o> -e <stderrfile.e> submit_manager.sh 
#echo $(hostname) > hostfile

. "/data/homezvol0/lilyw7/miniconda3/etc/profile.d/conda.sh"
conda activate psiresp-3.8


./01_test_energy.py  $MOLFILE