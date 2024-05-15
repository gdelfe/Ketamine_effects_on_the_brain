**Functions included in utilities_ketamine_analysis_v8**

**set_pub_plots** -\> set the output style of figures.

**bin_smooth_spk** -\> Given a spiking dataframe (df, with columns ['unit_id', 'spk’]), and the binsize (samples per second), it outputs a dataframe of spiking activity binned at the chosen time resolution (binsize) and after applying a convolution of chosen window.

**bin_spk** -\> Same as in **bin_smooth_spk**  but without applying a convolution that smooths the spiking activity.

**NormalizeData** -\> Function that normalizes a given input.

**kendall_tau_parallel** -\> Given a matrix input of cells x time, it computes (using parallelization) the Kendall tau correlation of each cell pair. It outputs three variables: a correlation matrix (cells x cells), a list of the cell pair IDS, and a list containing the kendall tau statistics and the corresponding p-value for each corresponding cell pair ID. Transposing the input results in the time-to-time correlations.

**pearson_parallel** -\> Given a matrix input of cells x time, it computes (using parallelization) the Pearson correlation of each cell pair. It outputs three variables: a correlation matrix (cells x cells), a list of the cell pair IDS, and a list containing the kendall tau statistics and the corresponding p-value for each corresponding cell pair ID. Transposing the input results in the time-to-time correlations.

**isomap_proj** -\> Given an input of behavior (timeseries), spiking activity (cells x time), the output dimensions (ndim), and the number of neighbors parameter (nneigh), this function computes the isomap projection and produces a behaviorally-relevant visualization. The output is the isomap projection of the shape ndim x time.

**compute_betti** -\> Given a 10-D isomap projection (or user-defined high-dimensional data), it computes the betti barcodes 0,1 and 2. This function includes options to downsample the data for computing betti 2 barcdode amongst others. The output comprises the three Betti barcodes to characterize persistent cohomology of the input.

**UMAP_proj** -\> Given an input of spiking activity (cells X time), and a timeseries of a relevant behavioral variable, this function computes the UMAP projection and a visualization of the result.

**kendall_visualization_PCo** -\> Given two inputs containing the Kendall tau correlations for the same cell pairs in two different conditions, it creates a correlation scatterplot visualization and outputs the Population Coordination value.

**kendall_visualization_PCo_doses** -\> Given four inputs containing the Kendall tau correlations for the same cell pairs in four different conditions (that is baseline, dose1, dose2 and dose3), it creates a correlation scatterplot visualization and outputs the Population Coordination between each condition and baseline (first array).

**tau_distribution_timeseries** -\> Given an input list containing Kendall tau correlations between the cell pairs in HPC-HPC, HPC-PFC, and PFC-PFC at different time-windows, and an input mask of which cell pair belongs to each of the three brain region comparison (HPC-HPC, HPC-PFC, and PFC-PFC), this function computes the proportions and ratios of positive and negative correlations across the different time-windows, per brain region of interest.

**tau_corr_dynamics** -\> Similar to **tau_distribution_timeseries**. Given an input list containing Kendall tau correlations between the cell pairs in HPC-HPC, HPC-PFC, and PFC-PFC during the 4 conditions (baseline, dose1, dose2, and dose3), and an input mask of which cell pair belongs to each of the three brain region comparison (HPC-HPC, HPC-PFC, and PFC-PFC), this function computes the proportions and ratios of positive and negative correlations across dosage, per brain region of interest.

**piechart_plot** -\> Visualization of the evolution of positively and negatively cell pair correlations across conditions.

**getMeanBinnedVector_2d** -\> It takes a matrix of spikes by time (or vice versa) as well as the desired number of temporal bins to be compressed to. The function splits the original matrix into \#bins and computes the average firing rate of each unit per time bin.

**neighbor_sum** -\> Given a rate map input, this function computes the spatial coherence, a metric we use to determine place cell properties.

**place_info_content** -\> Given the occupancy matrix, and the rate maps, it computes the place information context, a metric we use to determine place cell properties.

**nan_gaussian_filter** -\> Given a rate map input and the sigma value, it outputs a smoothed place map.

**conv_nan** -\> This function performs a convolution ignoring NaNs.

**fisherZ** -\> Given an r value, it computes the Fisher’s Z score.

**loadPhysiology** -\> Function that inputs parameters and outputs physiology .mat file and additional extracted and formatted information (metadata).

**runsTest** -\> Given a timeseries and a threshold s0, it computes the runs test for random alternation of the timeseries.

**place_cell_properties** -\> Given inputs of spiking activity (list of cells, inside each element, spiking timeseries), behavior data (2-D array with x and y coordinates timeseries), the number of bins in each direction (binsxy), and the sampling rate of the inputs (fsrate), it computes two measures of placeness, the coherence and the place information, to determine objectively (comparing to 100 or 1000 shuffles) if a cell is a place cell or not. The output comprises the following variables: occupancy map matrix, a list of the firing rate maps (length of \#cells), a list of place information for each cell, a list of coherence for each cell, and a Boolean array indicating whether a cell was categorized as place cell or non-place cell.

**normalizer** -\> Function that normalizes a given input time series.

**fr_map_ang** -\> Given inputs of spiking activity (list of length \#cells, each element containing the timeseries of activity), the behavior data (a 2-D array of x, y coordinates timeseries, the number of spatial bins, and the sampling rate (fsrate), this function output the linearized firing rate maps for each cell.

**fr_map** -\> Given inputs of spiking activity (list of length \#cells, each element containing the timeseries of activity), the behavior data (a 2-D array of x, y coordinates timeseries, the number of spatial bins in each x, y direction, and the sampling rate (fsrate), this function output the 2-D firing rate maps for each cell.


**upsample_speed** -\> Given the speed timeseries and the LFP timeseries, this function upsamples the speed (lower sampling rate, 100 Hz) to the sampling rate of the LFP signal (2500 Hz). It outputs the upsampled speed timeseries.

**create_speed_mask** -\> Given the upsampled speed timeseries, a time window length in samples, and a speed threshold, this function creates a Boolean array indicating whether the mean speed at each time window is above the defined threshold.

**create_sta_mask** -\> Given the upsampled speed timeseries, a time window length in samples, and a speed threshold, this function creates a Boolean array indicating whether the speed at all time-points at each time window is below the defined threshold.


**Functions for Markov process computation**

**compute_transition_matrix2** -\> Given a timeseries, it computes the transition probability matrix, P.

**multiply** -\> Function that multiplies two matrices A and B.

**matrix_power** -\> Function that calculates the power of a matrix input, A to a power, p.

**findProbability** -\> Function that computes the probability of reaching state F at time T after starting from S. The inputs to this function are as follows, **findProbability(G, N, F, S, T)**
- G, is transition matrix
- N, is number of states
- F, where to go
- S, from where
- T, how many samples
