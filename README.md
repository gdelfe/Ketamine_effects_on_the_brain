# Ketamine effects on hippocampus and prefrontal cortex

### Recordings

Electrophysiology and free behavior were recorded simultaneously while RS-ketamine doses (3, 10 and 30 mg/kg) were injected at minutes 30, 60 and 90 of a 2-hour session.

Data comprises a saline control session plus two RS-ketamine sessions (RSK), two days apart, for 5 mice.

Two neuropixels probes were used to capture simultaneous local field (LFP) and action potentials (AP) in hippocampus (HPC) and prefrontal cortex (PFC).

The action potentials were sampled at 30 kHz and LFP was sampled at 2500 Hz by 384 recording electrodes spaced 20 um apart depthwise along each probe.

| individual | date | drug | doses (mpk) | session index |
| --- | --- | --- | --- | --- |
| M015 | 2022-07-27 2022-08-01 2022-08-03 | SAL RSK RSK | 0, 0, 0 3, 10, 30 3, 10, 30 | 0 2 3 |
| M016 | 2022-07-28 2022-08-04 2022-08-06 | SAL RSK RSK | 0, 0, 0 3, 10, 30 3, 10, 30 | 1 4 5 |
| M017 | 2022-08-08 2022-08-12 | SAL RSK | 0, 0, 0 3, 10, 30 | 6 8 |
| M018 | 2022-08-11 2022-08-13 2022-08-15 | SAL RSK RSK | 0, 0, 0 3, 10, 30 3, 10, 30 | 7 9 10 |
| M023 | 2022-09-12 2022-09-14 2022-09-16 | SAL RSK RSK | 0, 0, 0 3, 10, 30 3, 10, 30 | 15 16 17 |

(M017 has a second RSK recording during which the rig went down and the recording was split into 2 sections; the session is presently ignored.)

The HPC probes passed through hippocampal strata oriens, pyramidale, radiatum, and lacunosum-moleculare before crossing the hippocampal fissure into the upper and lower dentate layers.

Each stratum has a characteristic current-source density profile, making it possible in principle to estimate electrode depth from local field electrophysiology (see Brankack, Stewart, Fox, Brain Research, 1993, Figure 2).

Analysis so far has shown that ketamine generally attenuates current source density power at all frequencies except a particular bnad of medium-high gamma. This acute effect increases with increasing dosage.

We are interested in how acute effects differ between the first ("pre-treatment") and second ("post-treatment") exposure. Some analyses suggest that the same effect level can be reached at a lower dosage in the post-treatment than in the pre-treatment session.

### Usage

Dependencies (mac environment):
1. sshfs (to mount `monk` filesystem)
2. python3
3. python-tk (`brew install python-tk`)

Installation:
1. Clone this repository to local filesystem
2. Run `make env` to build python environment

To use repository:
1. Run `make notebook` to start jupyter notebook server
2. Run `make mount` to mount remote filesystem
3. Run `make unmount` to unmount filesystem when done

### Repository structure

`gino` contains code from Gino,
1. converts raw LFP to CSD and computes running speed/LFP artifact masks (python)
    - accepts raw neuropixel files (.bin and .meta) as input and produces CSD and mask files (.mat) as output
    - each second is considered as a trial having movement speed 'none', 'low' or 'high'
    - the masks identify low-, high-speed movement and LFP artifacts
2. computes CSD power-spectral density (PSD) for each minute of each epoch throughout the session (matlab)
    - the epochs are 'B', 'L', 'M', 'H', each takes central 20 of 30 minutes
    - accepts CSD and mask files (.mat) as input and produces session PSD file as output (.mat)
3. compares PSD for each epoch between RSK sessions 1 and 2 (matlab)
4. analyzes spike fields (matlab)

`simon` contains code from Simon, which does spike sorting, cofiring analysis, manifold visualization, and other things.

`kfx` is space for new code built on top of Gino and Simon's framework.

`notebooks` has separate spaces with notebooks from each person.

`ref` contains metadata files, spreadsheets, etc. for common reference.

`data` contains raw or working data and is not tracked by version control.
