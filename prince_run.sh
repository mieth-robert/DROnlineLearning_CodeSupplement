#!/bin/bash

#SBATCH --time=10:00:00
#SBATCH --job-name=fullrun1_a
#SBATCH --mail-type=END
#SBATCH --mail-user=rm4796@nyu.edu

module purge
module load mosek/8.1.0.64
module load julia/0.6.3

srun julia ./run_dronlinelearning.jl ./cases/no_network.jl
srun julia ./run_dronlinelearning.jl ./cases/only_flow.jl
srun julia ./run_dronlinelearning.jl ./cases/only_voltage.jl
srun julia ./run_dronlinelearning.jl ./cases/all_constraints.jl
srun julia ./run_dronlinelearning.jl ./cases/fully_integrated.jl
