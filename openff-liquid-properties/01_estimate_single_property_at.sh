#!/bin/bash
#SBATCH -A lilyw7
#SBATCH -p free-gpu          ## partition/queue name
#SBATCH --nodes=1            ## (-N) number of nodes to use
#SBATCH --ntasks=1           ## (-n) number of tasks to launch
#SBATCH -t 200
#SBATCH --gres=gpu:V100:1 
#SBATCH --output=%x.log%A

echo $CUDA_VISIBLE_DEVICES

. "/data/homezvol0/lilyw7/miniconda3/etc/profile.d/conda.sh"

export ELF10=false

conda activate polymetrizer


conda list


./01_estimate_single_property.py $SLURM_ARRAY_TASK_ID