#!/bin/bash
#SBATCH --job-name=csd      # Job name
#SBATCH --output=job.out    # Output file
#SBATCH --error=job.out     # Error file
#SBATCH --time=2:00:00      # Time limit hrs:min:sec
#SBATCH --mem=16G           # Memory required per node
#SBATCH --nodes=1           # Number of nodes per SLURM job
#SBATCH --ntasks=1          # Number of MPI tasks per node
#SBATCH --cpus-per-task=16  # Number of CPU cores per MPI task
#SBATCH --partition=ccn     # Partition name
#SBATCH -C ib-rome          # 128 cores, 1 TB RAM, InfiniBand

# usage: sbatch --array=0-27 slurm_job.sh

# Set up environment for SLURM submission
module -q purge                       # purge current modules
module -q load openmpi                # load openmpi

ROOT="/mnt/home/larend/ceph/results/20240519/${SLURM_ARRAY_JOB_ID}"

mkdir -p "${ROOT}"
mpirun -n 1 python mpi_job.py \
    --output-path "${ROOT}" \
    --params-file="jobs.json"
