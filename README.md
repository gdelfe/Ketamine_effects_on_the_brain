# Ketamine effects on the brain, with focus on the Hippocampus and the Pre-Frontal Cortex 

This repo contains Python and Matlab code for the analysis of the Ketamine's effects on the brain a mouse model, work done within Andr&eacute; Fenton Lab, at NYU.

The mice receive 3 subsequent and increasing doses of Ketamine, one every 30 min, in a resting state dynamics (no specific task). The animals are free to move without restrictions in a 2D arena. </br>

The signals are simultaneously recorded with two Neuropixel probes, one located in the Hippocampus (HPC) and the other in the Pre-Frontal Cortex (PFC). Both Local Field Potentials (LFPs) signal and spikes (single units) are analyzed. </br>

The Matlab folder contains codes for: the pre-processing of the LFP; Power Spectral Density estimates (multi-taper methods); Spectrograms estimates; </br>
The Python folder contains codes for: behavioural variable extraction; LFP instantaneous phase extraction at given frequency (using Morlet Wavelet); Place cells identification; and more...
