#!/bin/bash



for mol in 18 210 229 251 268 614 643 776 785 1445; do
    sbatch --export=MOL=${mol} -J ${mol}_oe_charge 03_charge_single.sh
done