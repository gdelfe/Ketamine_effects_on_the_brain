import numpy as np
import pandas as pd
import pywt


def wavelet_transform(csd, freq_min, freq_max, freq_bins):
    """
    csd: single-channel CSD timeseries as pd.Series
    freq_min: freq of largest-scale wavelet, Hz
    freq_max: freq of smallest-scale wavelet, Hz
    freq_bins: apply transform at this many scales log-spaced from min to max
    """
    T = 1 / 1250
    # 'cmor1.5-1.0': complex Morlet wavelet, 1.0 Hz bandwidth 1.5 Hz
    wavelet = pywt.ContinuousWavelet('cmor1.5-1.0')
    freqs = np.logspace(np.log10(freq_min), np.log10(freq_max), freq_bins)
    scales = pywt.frequency2scale(wavelet, freqs * T)
    coef, _ = pywt.cwt(csd, scales, wavelet, sampling_period=T)
    return pd.DataFrame(coef.T, index=csd.index, columns=freqs)


def extract_phase(coef):
    """
    coef: output of wavelet transform as pd.DataFrame
    """
    phi = np.angle(coef.values)
    return pd.DataFrame(phi, index=coef.index, columns=coef.columns)


def resample_1250_250_Hz(df):
    """
    df: dataframe indexed at 1250 Hz in seconds
    """
    idx = df.index[::5]
    t1 = pd.date_range(start=0, periods=len(df), freq='0.8ms')
    t2 = pd.date_range(start=0, periods=len(df) // 5, freq='4ms')
    df = df.set_index(t1).reindex(t2).interpolate(method='linear')
    return df.set_index(idx)
