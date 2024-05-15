import h5py
import numpy as np
import pandas as pd
from scipy.io import loadmat


def load_lfp(path, startmin=0, endmin=1, offset=0, channels=None):
    """ Load raw voltages as recorded by the neuropixels probe. 

    path = '2022-08-11-01-55-00_M018_SAL_mPFC_HPC_0_0_0mpk_g0_t0.imec1.lf.bin'
    df = load_lfp(path)
    """
    
    with open(path.replace('.bin', '.meta'), 'r') as f:
        meta = {}
        for line in f.readlines():
            k, v = line.split('=')
            meta[k.strip('~')] = v.strip()

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
    int2volt = 0.6 / 512 / gain # convert each 2-byte integer to Volts
    if channels is None:
        # probe has 385 channels, the last one is not used for recording
        channels = np.arange(384)

    # open a memory map to the raw binary file and read out what we need
    ncols = nbytes // 2 // 385 - offset
    startbyte = offset * 2 * 385
    vals = np.memmap(path, dtype='int16', mode='r', shape=(385, ncols), offset=startbyte, order='F')
    start = int(startmin * 60 * fs)
    end = int(endmin * 60 * fs)
    voltage = vals[channels, start:end].T * int2volt  # [V]

    time = pd.Series(np.arange(start, end) / fs, name='time')  # [s]
    chan = pd.Series(channels, name='channel')
    return pd.DataFrame(voltage, index=time, columns=chan)


def load_csd(lfp, mask):
    """
    Load current source density with trial-wise speed information.

    path1 = 'lfp_epoch_all_trials_CSD.mat'
    path2 = 'mask_low_high_speed_CSD.mat'
    df = load_csd(path1, path2)
    """
    epochs = pd.Series(['B', 'L', 'M', 'H'], name='epoch')
    
    # len 4, each epoch is (20 minutes, 60 * 1250 ts, n sites)
    blob1 = loadmat(lfp)['lfp_all'][0][0]
    nsites = blob1[0].shape[-1]

    # len 8, first 4 are low-speed masks, last 4 are high-speed masks
    # 1: low/high speed, 0: not low/high speed or LFP artifact present
    blob2 = loadmat(mask)['mask'][0][0]

    # put LFP samples on the rows, channels on the columns
    minutes = pd.Series(range(20), name='minute')
    seconds = pd.Series(range(60), name='second')
    step = pd.Series(range(1250), name='step')
    chan = pd.Series(range(nsites), name='channel')
    idx = pd.MultiIndex.from_product([epochs, minutes, seconds, step])

    X = np.zeros((len(idx), nsites))
    l = []
    for i, epoch in enumerate(epochs):
        for minute in minutes:
            start = (i * 20 + minute) * 60 * 1250
            end = start + 60 * 1250
            X[start:end, :] = blob1[i][minute, :, :]

            # each second is assigned trial type 'none', 'low', 'high'
            trials = np.full(60, 'none')
            low = blob2[i][minute]
            hi = blob2[i + 4][minute]
            trials[np.all(low, axis=1)] = 'low'
            trials[np.all(hi, axis=1)] = 'high'
            for second in seconds:
                l.extend(1250 * [trials[second]])

    # index by tuple epoch, minute, second, trial, time
    df = pd.DataFrame(X, index=idx, columns=chan)
    idx = df.index.to_frame()
    idx.insert(3, 'trial', l)
    idx.insert(4, 'time', idx.second + idx.step / 1250)
    df.index = pd.MultiIndex.from_frame(idx.drop(columns='step'))
    return df


def load_psd(path):
    """
    Load power spectral density computed for each minute of each epoch.
    
    df = load_csd('psd_hpc_csd.mat')
    """
    blob = h5py.File(path)['psd']['HPC']
    epochs = pd.Series(['B', 'L', 'M', 'H'], name='epoch')
    minute = pd.Series(range(20), name='minute')
    freq = pd.Series(blob['f'][:, 0], name='frequency')
    idx = pd.MultiIndex.from_product([epochs, minute])
    df = pd.DataFrame(index=idx, columns=freq)
    for epoch in epochs:
        df.loc[epoch, :] = np.array(blob[epoch]).T
    return df