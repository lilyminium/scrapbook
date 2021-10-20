#!/bin/bash

mol=18
index=16
sbatch --export=MOL=${mol},ENTRY=${index} -J ${mol}_${index}_oe_eval 03_submit_single.sh

mol=251
index=255
sbatch --export=MOL=${mol},ENTRY=${index} -J ${mol}_${index}_oe_eval 03_submit_single.sh

mol=776
index=775
sbatch --export=MOL=${mol},ENTRY=${index} -J ${mol}_${index}_oe_eval 03_submit_single.sh



# for mol in 18 210 229 251 268 614 643 776 785 1445; do
#     sbatch --export=MOL=${mol} -J ${mol}_oe_eval 03_submit_single.sh
# done