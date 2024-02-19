#!/bin/bash -l
#SBATCH --job-name="VASP"
#SBATCH --account=s1167
#SBATCH --partition=normal
# #SBATCH --partition=debug
#SBATCH --time=24:00:00
#SBATCH --ntasks=2
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH --constraint=gpu

#SBATCH --output=log.out
#SBATCH --error=log.err


# runfile="$(printf 'run-%0.5i.sh' $id)"
# rundir="$(printf 'run-%0.5i' $id)"



echo "SLURM_NTASKS: $SLURM_NTASKS"
echo "SLURM_NTASKS_PER_NODE: $SLURM_NTASKS_PER_NODE"
echo "SLURM_CPUS_PER_TASK: $SLURM_CPUS_PER_TASK"

# load modules and run simulation
module load daint-gpu
module load QuantumESPRESSO
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export NO_STOP_MESSAGE=1
export CRAY_CUDA_MPS=1

export ESPRESSO_PSEUDO=$HOME/qe_pseudopotentials
# export ASE_ESPRESSO_COMMAND="srun pw.x"

idx=0

SCRIPT_DIR="/users/tayfurog/qe_interface/"
# PYTHON_DIR="/users/tayfurog/miniconda3/bin/"

python $SCRIPT_DIR/runQEArray.py\
	-calc_type freq_dfpt -geoms_path qe_opt_lowest_10_isolated_24atoms.extxyz -func PBE -idx $idx -memory 128
