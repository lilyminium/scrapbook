#!/bin/bash
#SBATCH -A lilyw7
#SBATCH -p free-gpu          ## partition/queue name
#SBATCH --nodes=1            ## (-N) number of nodes to use
#SBATCH --ntasks=1           ## (-n) number of tasks to launch
#SBATCH -t 100
#SBATCH --gres=gpu:V100:1 
#SBATCH --output=%x.log%A

echo $CUDA_VISIBLE_DEVICES

. "/data/homezvol0/lilyw7/miniconda3/etc/profile.d/conda.sh"

export OE_LICENSE=/data/homezvol0/lilyw7/oe_license.txt
export ELF10=true

conda activate polymetrizer


conda list


./01_estimate_single_property.py 1 #$SLURM_ARRAY_TASK_ID