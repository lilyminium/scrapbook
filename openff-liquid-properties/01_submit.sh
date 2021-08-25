#!/bin/bash

sbatch --array=1-4 01_estimate_single_property_oe.sh
sbatch --array=1-4 01_estimate_single_property_oe_elf10.sh
sbatch --array=1-4 01_estimate_single_property_at.sh

# sbatch --array=1-16925 01_estimate_single_property_oe.sh
# sbatch --array=1-16925 01_estimate_single_property_oe_elf10.sh
# sbatch --array=1-16925 01_estimate_single_property_at.sh