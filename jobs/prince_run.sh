#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --time=6:00:00
#SBATCH --mem=20GB
#SBATCH --job-name=2E6Soff
#SBATCH --mail-type=END
#SBATCH --mail-user=jipkim@nyu.edu

module purge
module load mosek/8.1.0.64
module load julia/0.6.3

srun julia run_dronlinelearning.jl cases/no_network.jl
srun julia run_dronlinelearning.jl cases/only_flow.jl
srun julia run_dronlinelearning.jl cases/only_voltage.jl
srun julia run_dronlinelearning.jl cases/all_constraints.jl
srun julia run_dronlinelearning.jl cases/fully_integrated.jl
