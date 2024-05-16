#!/bin/bash
#SBATCH --job-name=monkcopy    # Job name
#SBATCH --output=monkcopy.out  # Output file
#SBATCH --error=monkcopy.err   # Error file
#SBATCH --time=12:00:00        # Time limit hrs:min:sec
#SBATCH --partition=ccn        # Partition name
#SBATCH --nodes=1              # Number of nodes
#SBATCH --ntasks=1             # Number of tasks (typically 1 for a single command)
#SBATCH --cpus-per-task=16     # Number of CPU cores per task
#SBATCH --mem=16G              # Memory required per node

# monkcopy
# --------
# copy raw data files from monk to flatiron
# store creds in file `monkauth`
# link is ~100 MB/s
#
# usage: ./monkcopy.sh           (bash)
#        sbatch monkcopy.sh      (slurm)

# which kind of file to copy?
sig="ap" # 'ap' | 'lf'
ext="bin" # 'meta' | 'bin'

# these two raw files break naming convention of all others
sshpass -f monkauth scp -v "lukea@monk.cns.nyu.edu:/f/fentonlab/RAWDATA/NeuroPix/Ketamine/2022-07-27-07-41-00_M015_SAL_PFC_HPC_0_0_0mpk/07-27-2022-M015_SAL_0.09_0.3_0.9_g0_imec1/07-27-2022-M015_SAL_0.09_0.3_0.9_g0_t0.imec1.${sig}.${ext}" "2022-07-27-07-41-00_M015_SAL_PFC_HPC_0_0_0mpk_g0_t0.imec1.${sig}.${ext}"

sshpass -f monkauth scp -v "lukea@monk.cns.nyu.edu:/f/fentonlab/RAWDATA/NeuroPix/Ketamine/2022-07-28-13-19-00_M016_SAL_PFC_HPC_0_0_0mpk/07-28-2022_M016_SAL_0.09_0.3_0.9_g0_imec1/07-28-2022_M016_SAL_0.09_0.3_0.9_g0_t0.imec1.${sig}.${ext}" "2022-07-28-13-19-00_M016_SAL_PFC_HPC_0_0_0mpk_g0_t0.imec1.${sig}.${ext}"

# these three raw files have an underbar instead of hyphen in one spot
sshpass -f monkauth scp -v "lukea@monk.cns.nyu.edu:/f/fentonlab/RAWDATA/NeuroPix/Ketamine/2022-08-01_04-30-00_M015_RSK_mPFC_HPC_3_10_30mpk/2022-08-01_04-30-00_M015_RSK_mPFC_HPC_3_10_30mpk_g0_imec1/2022-08-01_04-30-00_M015_RSK_mPFC_HPC_3_10_30mpk_g0_t0.imec1.${sig}.${ext}" "2022-08-01-04-30-00_M015_RSK_mPFC_HPC_3_10_30mpk_g0_t0.imec1.${sig}.${ext}"

sshpass -f monkauth scp -v "lukea@monk.cns.nyu.edu:/f/fentonlab/RAWDATA/NeuroPix/Ketamine/2022-08-03_12-45-00_M015_RSK_mPFC_HPC_3_10_30mpk/2022-08-03_12-45-00_M015_RSK_mPFC_HPC_3_10_30mpk_g0_imec1/2022-08-03_12-45-00_M015_RSK_mPFC_HPC_3_10_30mpk_g0_t0.imec1.${sig}.${ext}" "2022-08-03-12-45-00_M015_RSK_mPFC_HPC_3_10_30mpk_g0_t0.imec1.${sig}.${ext}"

sshpass -f monkauth scp -v "lukea@monk.cns.nyu.edu:/f/fentonlab/RAWDATA/NeuroPix/Ketamine/2022-08-04_03-45-00_M016_RSK_mPFC_HPC_3_10_30mpk/2022-08-04_03-45-00_M016_RSK_mPFC_HPC_3_10_30mpk_g0_imec1/2022-08-04_03-45-00_M016_RSK_mPFC_HPC_3_10_30mpk_g0_t0.imec1.${sig}.${ext}" "2022-08-04-03-45-00_M016_RSK_mPFC_HPC_3_10_30mpk_g0_t0.imec1.${sig}.${ext}"

# this session was broken in two since the behavior rig went down
if [ "$ext" = "meta" ]; then

    # metadata for first chunk
    sshpass -f monkauth scp -v "lukea@monk.cns.nyu.edu:/f/fentonlab/RAWDATA/NeuroPix/Ketamine/2022-08-10-03-33-00_M017_RSK_mPFC_HPC_3_10_30mpk/2022-08-10-03-33-00_M017_RSK_mPFC_HPC_3_10_30mpk_g0_imec1/2022-08-10-03-33-00_M017_RSK_mPFC_HPC_3_10_30mpk_g0_t0.imec1.${sig}.${ext}" "2022-08-10-03-33-00_M017_RSK_mPFC_HPC_3_10_30mpk_g0_t0.imec1.${sig}1.${ext}"

    # metadata for second chunk
    sshpass -f monkauth scp -v "lukea@monk.cns.nyu.edu:/f/fentonlab/RAWDATA/NeuroPix/Ketamine/2022-08-10-03-33-00_M017_RSK_mPFC_HPC_3_10_30mpk_2/2022-08-10-03-33-00_M017_RSK_mPFC_HPC_3_10_30mpk_2_g0_imec1/2022-08-10-03-33-00_M017_RSK_mPFC_HPC_3_10_30mpk_2_g0_t0.imec1.${sig}.${ext}" "2022-08-10-03-33-00_M017_RSK_mPFC_HPC_3_10_30mpk_g0_t0.imec1.${sig}2.${ext}"

    # metadata for combined is just a copy of first chunk metadata
    sshpass -f monkauth scp -v "lukea@monk.cns.nyu.edu:/f/fentonlab/RAWDATA/NeuroPix/Ketamine/2022-08-10-03-33-00_M017_RSK_mPFC_HPC_3_10_30mpk_all/2022-08-10-03-33-00_M017_RSK_mPFC_HPC_3_10_30mpk_g0_imec1/2022-08-10-03-33-00_M017_RSK_mPFC_HPC_3_10_30mpk_g0_t0.imec1.${sig}.${ext}" "2022-08-10-03-33-00_M017_RSK_mPFC_HPC_3_10_30mpk_g0_t0.imec1.${sig}.${ext}"

# else

    # first chunk of recording
    sshpass -f monkauth scp -v "lukea@monk.cns.nyu.edu:/f/fentonlab/RAWDATA/NeuroPix/Ketamine/2022-08-10-03-33-00_M017_RSK_mPFC_HPC_3_10_30mpk/2022-08-10-03-33-00_M017_RSK_mPFC_HPC_3_10_30mpk_g0_imec1/2022-08-10-03-33-00_M017_RSK_mPFC_HPC_3_10_30mpk_g0_t0.imec1.${sig}.${ext}" "2022-08-10-03-33-00_M017_RSK_mPFC_HPC_3_10_30mpk_g0_t0.imec1.${sig}1.${ext}"

    # second chunk of recording
    sshpass -f monkauth scp -v "lukea@monk.cns.nyu.edu:/f/fentonlab/RAWDATA/NeuroPix/Ketamine/2022-08-10-03-33-00_M017_RSK_mPFC_HPC_3_10_30mpk_2/2022-08-10-03-33-00_M017_RSK_mPFC_HPC_3_10_30mpk_2_g0_imec1/2022-08-10-03-33-00_M017_RSK_mPFC_HPC_3_10_30mpk_2_g0_t0.imec1.${sig}.${ext}" "2022-08-10-03-33-00_M017_RSK_mPFC_HPC_3_10_30mpk_g0_t0.imec1.${sig}2.${ext}"

    # first and second chunks concated (LFP didn't stop recording, just behavior)
    sshpass -f monkauth scp -v "lukea@monk.cns.nyu.edu:/f/fentonlab/RAWDATA/NeuroPix/Ketamine/2022-08-10-03-33-00_M017_RSK_mPFC_HPC_3_10_30mpk_all/2022-08-10-03-33-00_M017_RSK_mPFC_HPC_3_10_30mpk_g0_imec1/2022-08-10-03-33-00_M017_RSK_mPFC_HPC_3_10_30mpk_g0_t0.imec1.${sig}2.${ext}" "2022-08-10-03-33-00_M017_RSK_mPFC_HPC_3_10_30mpk_g0_t0.imec1.${sig}.${ext}"

fi

# do the rest
for sess in \
2022-08-13-03-57-00_M018_RSK_mPFC_HPC_3_10_30mpk \
2022-08-15-02-25-00_M018_RSK_mPFC_HPC_3_10_30mpk \
2022-08-06-05-20-00_M016_RSK_mPFC_HPC_3_10_30mpk \
2022-08-08-04-05-00_M017_SAL_mPFC_HPC_0_0_0mpk \
2022-08-11-01-55-00_M018_SAL_mPFC_HPC_0_0_0mpk \
2022-08-12-02-35-00_M017_RSK_mPFC_HPC_3_10_30mpk \
2022-09-12-04-40-00_M023_SAL_mPFC_HPC_0_0_0mpk \
2022-09-14-03-35-00_M023_RSK_mPFC_HPC_3_10_30mpk \
2022-09-16-02-10-00_M023_RSK_mPFC_HPC_3_10_30mpk
do
    sshpass -f monkauth scp -v "lukea@monk.cns.nyu.edu:/f/fentonlab/RAWDATA/NeuroPix/Ketamine/${sess}/${sess}_g0_imec1/${sess}_g0_t0.imec1.${sig}.${ext}" "${sess}_g0_t0.imec1.${sig}.${ext}"
done
