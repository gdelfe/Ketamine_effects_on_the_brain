import numpy as np
import pandas as pd
import quantities as pq
from scipy.signal import butter, iirnotch, filtfilt, lfilter

from kfx.icsd import SplineiCSD


def clean_lfp(lfp):
    # 374, 375, 376, 377, 378, 379, 380, 381, 382, 383 don't carry signal
    chan = np.setdiff1d(lfp.columns, [383, 376, 380, 378, 379, 382, 377, 381, 374, 375])
    # 191 is a dead channel (replace it with 190)
    chan = np.where(chan == 191, 190, chan)
    return lfp[chan]


def preprocess_lfp(lfp):
    X = lfp.values

    # subtract DC offset for each channel
    X = X - X.mean(axis=0)
    X = (X[:, :-1:2] + X[:, 1::2]) / 2
    # bandpass filter between 1 and 300 Hz, notch filter out 60 Hz ground
    b, a = butter(N=5, Wn=[1 / 1250, 300 / 1250], btype='band')
    X = filtfilt(b, a, X, axis=0)
    b, a = iirnotch(w0=60, Q=200, fs=2500)
    X = lfilter(b, a, X, axis=0)
    # decimate to 1250 Hz
    X = X[::2]

    time = lfp.index[::2]
    site = pd.Series(lfp.columns[::2] // 2, name='site')
    return pd.DataFrame(X, index=time, columns=site)


def compute_csd(lfp):
    h = 20 * 1e-6 * pq.m            # inter-electrode spacing [m]
    d = 0.5 * 1e-3 * pq.m           # electrode diameter [m]
    sigma = 0.3 * pq.S / pq.m       # intracellular medium conductivity [S/m]
    sigma_top = 0.3 * pq.S / pq.m   # cortical surface conductivity [S/m]
    z = np.arange(lfp.shape[1]) * h # electrode coordinates [m]
    X = lfp.T.values * pq.V         # recorded potential [V]
    estimator = SplineiCSD(
        X, coord_electrode=z, diam=d, sigma=sigma, sigma_top=sigma_top,
        tol=1e-12, f_type='gaussian', f_order=(20, 5), num_steps=len(z)
    )
    csd = estimator.get_csd()
    csd = estimator.filter_csd(csd) # apply spatial smoothing
    return pd.DataFrame(csd.T, index=lfp.index, columns=lfp.columns)