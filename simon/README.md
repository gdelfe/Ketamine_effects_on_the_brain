**DOCUMENTATION FOR ANALYSIS CODE FOR HEAD-FIXED RECORDINGS UNDER KETAMINE**

This repository contains the analysis codes used to analyze neural activity and behavior at the single and population levels. The following analysis are included:

1) **Preprocessing of neural and behavioral data**
2) **Preliminary analysis**
3) **Single cell analysis**
4) **Population analysis (Cofiring and Population coordination)**
5) **Manifold analysis**
6) **Brain region identification (HPC)**
7) **Preprocessing of neural and behavioral data for simultaneous LFP and spiking analysis**


**Raw data nomenclature**

Filename: date + mouse ID + drug + brain region 1 + brain region 2 + dosage

Mouse ID (M015-M033)

Drug (SAL=saline, RSK=RS-ketamine, RK = R-ketamine, SK=S-ketamine, PCP=Phencyclidine, 1020, CP101606,

Example: 2022-08-06-05-20-00_M016_RSK_mPFC_HPC_3_10_30mpk_HPCprobe

This recording corresponds to mouse M16 recorded on August 6th, 2022, under a three-dosage protocol of 3, 10 and 30 mg/kg of RS-ketamine, injected in min 30, 60 and 90 of recording.

Behavior and synchronization files do not have a final handle. Electrophysiology files will have either a mPFC_probe or HPC_probe handle depending on the location in the brain of the probe.

Each recording produces the following set of raw data files (example below):

-   AP binary file (one per probe):

    2022-08-06-05-20-00_M016_RSK_mPFC_HPC_3_10_30mpk_g0_t0.ap.bin

    2022-08-06-05-20-00_M016_RSK_mPFC_HPC_3_10_30mpk_g0_t0.ap.bin

-   LF binary file (one per probe):

2022-08-06-05-20-00_M016_RSK_mPFC_HPC_3_10_30mpk_g1_t0.lf.bin

2022-08-06-05-20-00_M016_RSK_mPFC_HPC_3_10_30mpk_g1_t0.lf.bin

-   Spike sorted folder (one per probe):

2022-08-06-05-20-00_M016_RSK_mPFC_HPC_3_10_30mpk_HPCprobe

2022-08-06-05-20-00_M016_RSK_mPFC_HPC_3_10_30mpk_PFCprobe

-   Behavior file:

    2022-08-06-05-20-00_M016_RSK_mPFC_HPC_3_10_30mpk.mat

-   Sync file:

    2022-08-06-05-20-00_M016_RSK_mPFC_HPC_3_10_30mpk.mat

