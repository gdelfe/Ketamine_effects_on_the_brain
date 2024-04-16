# Ketamine effects on the brain, with focus on the Hippocampus and the Pre-Frontal Cortex 

This repo contains Python and Matlab code for the analysis of the Ketamine's effects on the brain a mouse model, work done within <a href= "https://www.fentonlab.com" > Andr&eacute; Fenton Lab </a>, at NYU.

The mice receive 3 subsequent and increasing doses of Ketamine, one every 30 min, in a resting state dynamics (no specific task). The animals are free to move without restrictions in a 2D arena. </br>

The signals are simultaneously recorded with two Neuropixel probes, one located in the Hippocampus (HPC) and the other in the Pre-Frontal Cortex (PFC). Both Local Field Potentials (LFPs) signal and spikes (single units) are analyzed. </br>

---

Dependencies (mac environment):
1. sshfs
2. python3
3. brew install python-tk

One time setup:
1. Clone this repository to local filesystem
2. Run 'make env' to build python environment

To use repository:
1. Run 'make mount' to mount remote filesystem
2. Run 'make notebook' to start jupyter notebook server
3. Run 'make unmount' to unmount filesystem when done

To set up local analysis:
1. Open `recid_CA1_DG_id_modified.xlsx` in Excel/Pages and export as `CSV recid_CA1_DG_id_modified.xlsx`
2. Change line 107 in `utils_general.py` to local path e.g. `in_file = '/Users/lukearend/phd/kfx/ref/recid_CA1_DG_id_modified.csv'`
3. Use `1-convert-paths.ipynb` to convert paths in `HPC_lfp_paths.file` and `PFC_lfp_paths.file` to paths on local filesystem, store local paths in files `HPC_lfp_paths_local.pkl` and `PFC_lfp_paths_local.pkl`
4. Change `Ketamine_LFP_CSD_preprocessing.ipynb` line 50-51 to point to local path files, e.g. `HPC_path_file = '/Users/lukearend/phd/kfx/ref/HPC_lfp_paths_local.pkl; PFC_path_file ='/Users/lukearend/phd/kfx/ref/PFC_lfp_paths_local.pkl'`
5. Copy directory `C://Users/fentonlab/Desktop/Gino/behaviour` to local path `data/behaviour`
6. Change `utils_general.py` line 140 to point to `data/behavior` e.g. `speed_path = '/Users/lukearend/phd/kfx/data/behaviour'`
