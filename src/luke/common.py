import os
import pandas as pd

sessions = pd.read_csv('/Users/lukearend/phd/kfx/ref/sessions.csv', index_col=0)
recordings = pd.read_csv('/Users/lukearend/phd/kfx/ref/recordings.csv', index_col=0)
metadata = pd.read_csv('/Users/lukearend/phd/kfx/ref/metadata.csv', index_col=[0, 1, 2])


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
