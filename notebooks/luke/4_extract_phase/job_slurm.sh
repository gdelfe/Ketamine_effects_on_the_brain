#!/bin/bash
#SBATCH --job-name=lukejob  # Job name
#SBATCH --output=job.out    # Output file
#SBATCH --error=job.out     # Error file
#SBATCH --time=1:00:00      # Time limit hrs:min:sec
#SBATCH --mem=128G          # Memory required per node
#SBATCH --nodes=1           # Number of nodes per SLURM job
#SBATCH --ntasks=8          # Number of MPI tasks per node
#SBATCH --cpus-per-task=16  # Number of CPU cores per MPI task
#SBATCH --partition=ccn     # Partition name
#SBATCH -C ib-rome          # 128 cores, 1 TB RAM, InfiniBand

module purge
module load openmpi

RESULTS_ROOT=$1
outpath="${RESULTS_ROOT}/${SLURM_ARRAY_JOB_ID}"

mkdir -p ${outpath}
mpirun -n 8 python job_mpi.py \
    --output-path ${outpath} \
    --params-file="params.json"
