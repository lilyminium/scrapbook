#!/bin/bash
#SBATCH -A lilyw7
#SBATCH -p free          ## partition/queue name
#SBATCH --nodes=1            ## (-N) number of nodes to use
#SBATCH --ntasks=1           ## (-n) number of tasks to launch
#SBATCH -t 100
#SBATCH --output=%x.log%A

echo $CUDA_VISIBLE_DEVICES

. "/data/homezvol0/lilyw7/miniconda3/etc/profile.d/conda.sh"
export OE_LICENSE=/data/homezvol0/lilyw7/oe_license.txt
export ELF10=false

conda activate psiresp-3.8


conda list


./03_generate_charges.py $MOL
