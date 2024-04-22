import numpy as np
import pandas as pd


def load_potential(path, startmin=0, endmin=1, offset=0, channels=None):
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
