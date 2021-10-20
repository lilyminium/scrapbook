#!/bin/bash
#SBATCH -A lilyw7
#SBATCH -p free-gpu          ## partition/queue name
#SBATCH --nodes=1            ## (-N) number of nodes to use
#SBATCH --ntasks=12           ## (-n) number of tasks to launch
#SBATCH -t 2000
#SBATCH --gres=gpu:V100:1 
#SBATCH --output=%x.log-%A

echo $CUDA_VISIBLE_DEVICES

. "/data/homezvol0/lilyw7/miniconda3/etc/profile.d/conda.sh"
export OE_LICENSE=/data/homezvol0/lilyw7/oe_license.txt
export ELF10=false

conda activate psiresp-3.8


conda list

DELTA=1.0
./05_write_input.py $DELTA
cd "05_vdw/${DELTA}"

# ./05_write_am1_input.py
# cd 05_vdw/openeye

echo "STARTING NONBONDED RUN ${DELTA}"
nonbonded optimization run --restart true
echo "STARTING NONBONDED ANALYZE ${DELTA}"
nonbonded optimization analyze