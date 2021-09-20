#!/bin/bash
#SBATCH -p scavenge
#SBATCH --qos=scavenger
#SBATCH --account=mobley
#SBATCH -n 1
#SBATCH -t 4000
#SBATCH --output=%x_%A.%a.log


. "/export/nfs0home/lilyw7/miniconda3/etc/profile.d/conda.sh"

export OE_LICENSE=/export/nfs0home/lilyw7/oe_license.txt

conda activate polymetrizer

conda list



./04c_generate_elf10.py $SLURM_ARRAY_TASK_ID