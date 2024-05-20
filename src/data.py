import h5py
import numpy as np
import pandas as pd
from scipy.io import loadmat
import src


def load_metadata(path):
    """ Load the metadata associated with a neuropixels binary. """
    path = str(path).replace('.bin', '.meta')
    meta = {}
    with open(path, 'r') as f:
        for line in f.readlines():
            k, v = line.split('=')
            meta[k.strip('~')] = v.strip()
    return meta


def load_ephys(path, startmin=0, endmin=1, offset=0, channels=None, cleaned=False):
    """ Load raw voltag  recorded by the neuropixels probe.

    path = '2022-08-11-01-55-00_M018_SAL_mPFC_HPC_0_0_0mpk_g0_t0.imec1.lf.bin'
    df = load_lfp(path)
    """
    meta = load_metadata(path)

    # sanity check the input file
    assert int(meta['nSavedChans']) == 385
    assert int(meta['imSampRate']) in [2500, 30000]
    assert int(meta['imMaxInt']) == 512
    assert float(meta['imAiRangeMax']) == 0.6
    assert float(meta['imAiRangeMin']) == -0.6
    gain = meta['imroTbl'].split(sep=')')[1:-1]
    gain = np.array([int(i.split(sep=' ')[4]) for i in gain])
    assert len(np.unique(gain)) == 1
    gain = gain[0]

    nbytes = int(meta['fileSizeBytes'])
    fs = int(meta['imSampRate'])
    int2volt = 0.6 / 512 / gain # convert each 2-byte integer to volts
    if channels is None:
        # probe has 385 channels, the last one is not used for recording
        channels = np.arange(384)

    # open a memory map to the raw binary file and read out what we need
    ncols = nbytes // 2 // 385 - offset
    startbyte = offset * 2 * 385
    vals = np.memmap(path, dtype='int16', mode='r', shape=(385, ncols), offset=startbyte, order='F')
    start = int(startmin * 60 * fs)
    end = int(endmin * 60 * fs)
    vals = vals[channels, start:end].T
    voltage = vals * int2volt # [V]
    
    # in all HPC and PFC recordings, channel 191 is dead
    if cleaned:
        voltage[:, 191] = voltage[:, 190]

    time = pd.Series(np.arange(start, end) / fs, name='time')  # [s]
    chan = pd.Series(channels, name='channel')
    return pd.DataFrame(voltage, index=time, columns=chan)


def load_behavior(recnum):
    datadir = src.paths.DATA / 'behavior'
    
    x = np.load(datadir / 'x_aln.npy', allow_pickle=True)[recnum]
    y = np.load(datadir / 'y_aln.npy', allow_pickle=True)[recnum]
    yaw = np.load(datadir / 'yaw_aln.npy', allow_pickle=True)[recnum]
    speed = np.load(datadir / 'speed_aln.npy', allow_pickle=True)[recnum]
    behavior = {'x': x, 'y': y, 'yaw': yaw, 'speed': speed}
    
    time = np.arange(len(x)) / 100
    return pd.DataFrame(behavior, index=time)
