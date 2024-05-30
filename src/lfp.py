import numpy as np
import pandas as pd
import quantities as pq
from scipy.signal import butter, iirnotch, filtfilt, lfilter
from src.icsd import SplineiCSD


def remove_dc_offset(X):
    """ Remove DC offset channel-wise """
    return X - X.mean(axis=0)


def combine_neighbors(X):
    """
    Applied once, this averages over site pairs at the same depth.
    After that, it averages adjacent sites depthwise along the probe.
    """
    return (X[:, :-1:2] + X[:, 1::2]) / 2


def bandpass_1_300Hz(X):
    nyq = 2500 / 2
    b, a = butter(N=5, Wn=[1 / nyq, 300 / nyq], btype='band')
    return filtfilt(b, a, X, axis=0)


def notch_60Hz(X):
    """ Filter out ground noise """
    b, a = iirnotch(w0=60, Q=200, fs=2500) # don't change Q
    return lfilter(b, a, X, axis=0)


def decimate_by_2(X):
    return X[::2]


def compute_csd(X):
    """ Estimate current source density from LFP. """
    h = 20 * 1e-6 * pq.m             # inter-electrode spacing [m]
    d = 0.5 * 1e-3 * pq.m            # electrode diameter [m]
    sigma = 0.3 * pq.S / pq.m        # intracellular medium conductivity [S/m]
    sigma_top = 0.3 * pq.S / pq.m    # cortical surface conductivity [S/m]
    z = np.arange(X.shape[1]) * h    # electrode coordinates [m]
    X = X.T * pq.V                   # recorded potential [V]
    estimator = SplineiCSD(
        X, coord_electrode=z, diam=d, sigma=sigma, sigma_top=sigma_top,
        tol=1e-12, f_type='gaussian', f_order=(20, 5), num_steps=len(z)
    )
    csd = estimator.get_csd()
    csd = estimator.filter_csd(csd)  # apply spatial smoothing
    return csd.T
