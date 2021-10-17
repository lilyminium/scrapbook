#!/bin/bash

for molfile in $(ls input/*.json) ; do
    sbatch --export=MOLFILE=$molfile -J $molfile 01_hf_run_single.sh
done


