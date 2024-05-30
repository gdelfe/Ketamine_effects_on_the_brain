#!/bin/bash
source ~/kfx/env/bin/activate
sbatch --array=0-335 job_slurm.sh ~/ceph/results/20240530/
