# Ketamine effects on the brain, with focus on the Hippocampus and the Pre-Frontal Cortex 

This repo contains Python and Matlab code for the analysis of the Ketamine's effects on the brain a mouse model, work done within <a href= "https://www.fentonlab.com" > Andr&eacute; Fenton Lab </a>, at NYU.

The mice receive 3 subsequent and increasing doses of Ketamine, one every 30 min, in a resting state dynamics (no specific task). The animals are free to move without restrictions in a 2D arena. </br>

The signals are simultaneously recorded with two Neuropixel probes, one located in the Hippocampus (HPC) and the other in the Pre-Frontal Cortex (PFC). Both Local Field Potentials (LFPs) signal and spikes (single units) are analyzed. </br>

---

One time setup:
1. Clone this repository to local filesystem
2. Run 'make env' to build python environment

To use repository:
1. Run 'make mount' to mount remote filesystem
2. Run 'make notebook' to start jupyter notebook server
3. Run 'make unmount' to unmount filesystem when done
