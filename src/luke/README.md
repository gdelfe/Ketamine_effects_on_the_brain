To set up local analysis:
1. Open `recid_CA1_DG_id_modified.xlsx` in Excel/Pages and export as `CSV recid_CA1_DG_id_modified.xlsx`
2. Change line 107 in `utils_general.py` to local path e.g. `in_file = '/Users/lukearend/phd/kfx/ref/recid_CA1_DG_id_modified.csv'`
3. Use `1-comiple-paths.ipynb` to compile paths in `HPC_lfp_paths.file` and `PFC_lfp_paths.file` to a CSV file containing all session metadata.
4. Change `Ketamine_LFP_CSD_preprocessing.ipynb` line 50-51 to point to local path files, e.g. `HPC_path_file = '/Users/lukearend/phd/kfx/ref/HPC_lfp_paths_local.pkl; PFC_path_file ='/Users/lukearend/phd/kfx/ref/PFC_lfp_paths_local.pkl'`
5. Copy directory `C://Users/fentonlab/Desktop/Gino/behaviour` to local path `data/behaviour`
6. Change `utils_general.py` line 140 to point to `data/behavior` e.g. `speed_path = '/Users/lukearend/phd/kfx/data/behaviour'`
