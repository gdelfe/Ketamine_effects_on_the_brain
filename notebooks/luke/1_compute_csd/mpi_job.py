import argparse
import json
import logging
import os
import sys
from mpi4py import MPI
from job import main as run_job


def main(output_path, params_file):
    sjobid = int(os.environ['SLURM_ARRAY_JOB_ID'])
    taskid = int(os.environ['SLURM_ARRAY_TASK_ID'])
    ntasks = int(os.environ['SLURM_ARRAY_TASK_COUNT'])
    
    comm = MPI.COMM_WORLD
    threadid = comm.Get_rank()
    nthreads = comm.Get_size()
    
    logging.basicConfig(
        level=logging.INFO,
        format=f'[%(asctime)s][job {sjobid}][task {taskid}][thread {threadid}] %(message)s',
        handlers=[logging.StreamHandler()]
    )
    
    idx = threadid + (nthreads * taskid)
    logging.info(f'received job index {idx}')
    
    with open(params_file) as f:
        jobs = json.load(f)
    if idx > len(jobs) - 1:
        logging.info(f'job index {idx} is out of range, exiting')
        sys.exit(0)
    
    params = jobs[idx]
    run_job(output_path, params)

    
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--output-path', required=True)
    parser.add_argument('-p', '--params-file', required=True)
    args = parser.parse_args()
    main(args.output_path, args.params_file)
