# **Brain region identification (HPC)**

An estimate of the probe channel at which the different hippocampal subregions are located is computed using the function **HPC_lfp_CA1_DG_fisure_detection_v3**. The following computations are performed inside this function:

1.  Obtain channels with high theta power (**get_theta_ch_v2**).
2.  Filter LFP signal in the theta band (**theta_fir).**
3.  Apply spline smoothing to the theta band LFP signal to find channels where the theta phase reverses (i.e., location of the fissure).
4.  Plot Current Source Density filtered LFP together with smoothed theta band LFP signal (**plot_lfps_theta, theta_fir, csd_fir**)
5.  Obtain channels with high ripple power (**get_swr_ch_v2**).
6.  Print summary of estimated subregion locations (theta band and ripple channels= CA1, phase reversal= fissure, and Dentate Girus below).

Additional functions used during this computation:

**bandpass_filter** -\> Given a signal, it applies a bandpass filter between the user-defined low and high cut frequencies, and outputs the filtered signal.

**notch_filter** -\> Given a signal and a notch frequency, this function applies a band-rejection filter (notch filter), and outputs the filtered signal.
