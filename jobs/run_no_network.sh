#!/bin/bash

#SBATCH --job-name=fullrun1
#SBATCH --mail-type=END
#SBATCH --mail-user=rm4796@nyu.edu
#SBATCH --output=no_network.out

module purge
module load mosek/8.1.0.64
module load julia/0.6.3

cd ..

srun julia ./run_dronlinelearning.jl ./cases/no_network.jl
