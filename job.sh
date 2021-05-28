#!/bin/bash
#SBATCH -n 16
module load gromacs/5.1.4

goo-job-nanny srun mdrun -v -s md_0.tpr -o md_0.trr -e md_0.edr -c md_0.gro -g md_0.log
