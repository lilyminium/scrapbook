#!/bin/bash

for molfile in $(ls input/*.json) ; do
    for state in gas water ; do
        sbatch --export=MOLFILE=${molfile},STATE=${state} -J ${state}_${molfile} 01_resp2_single.sh
    done
done