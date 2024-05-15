
Code is included in Ketamine_HPC_spkforLFP_preprocessing_v6 Jupyter notebook.

# **Preprocessing of neural and behavioral data for simultaneous LFP and spiking analysis**

1) Search for all the filenames in the directories containing raw data: behavior, sync, ephys files.
2) Search for LFP raw files.
3) Load behavior and sync files. Extract alignment information (timestamps for ephys and timestamps for behavior trimming)
4) Load spike sorted files for both HPC and PFC probes, extract spiking activity at a chosen time resolution (in this case, finer resolution, 250 Hz), and trimmed according to alignment timestamps
5) Bin aligned behavior variables to the same time resolution as spiking activity.
6) Save for analysis of Freqency-Phase Place Cell analysis.
