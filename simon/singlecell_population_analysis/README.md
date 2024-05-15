
Code is included in Ketamine_singlecell_population_analysis Jupyter notebook.

# **Preprocessing of neural and behavioral data**

1) Search for all the filenames in the directories containing raw data: behavior, sync, ephys files.
2) Load behavior and sync files. Extract alignment information (timestamps for ephys and timestamps for behavior trimming)
3) Load spike sorted files for both HPC and PFC probes, extract spiking activity at a chosen time resolution and trimmed according to alignment timestamps
4) Bin aligned behavior variables to the same time resolution as spiking activity

# **Single cell analysis**

5) Behavior analysis: 2-D trajectories, speed and moving vs stationary masking
6) Preliminary firing rate analysis per condition
7) Place cell properties in both HPC and PFC single units. There are many choices to do here: per condition (dose), within conditions, moving vs stationary, etc.
8) Linearized firing rate maps across conditions
9) Markov process - Return time to spatial bin

# **Population analysis (Cofiring and Population coordination)**

10) Concatenate recordings (cell-wise) from both HPC and PFC.
11) Compute pairwise Kendall Tau correlations for all cell pairs, that is HPC-HPC pairs, PFC-PFC pairs and HPC-PFC pairs. This can be found for a number of conditions such as all timepoints, only moving, only stationary, and final comparison using Population coordination across the tested conditions.
12) Compute AV (Activity Vectors) correlation matrices

# **Manifold analysis**

13) Fit Isomap projections and UMAP projections.
