import json
import logging
import time
import numpy as np
import pandas as pd
import src.paths
from src.data import load_ephys, load_behavior
from src.lfp import remove_dc_offset, combine_neighbors, decimate_by_2
from src.lfp import bandpass_1_300Hz, notch_60Hz, compute_csd


def raw_lfp_to_csd(X):
    X = remove_dc_offset(X)
    X = combine_neighbors(X)
    X = bandpass_1_300Hz(X)
    X = notch_60Hz(X)
    X = decimate_by_2(X)
    X = compute_csd(X)
    X = combine_neighbors(X)
    return X


def compute_session_csd(date, region):
    with open(src.paths.REF / 'date_recoffset.json', 'r') as f:
        offsets = json.load(f)
    
    suffix = {'hpc': '.imec1.lf.bin', 'pfc': '.imec0.lf.bin'}[region]
    file = str(src.paths.session_path(date)) + suffix
    path = src.paths.DATA / date / file
    lfpstart, lfpend = offsets[date]
    duration = (lfpend - lfpstart) / 2500 / 60
    
    logging.info(f'path: {path}')
    logging.info(f'recording start, end: {lfpstart}, {lfpend}')
    logging.info(f'duration: {duration} min')

    intervals = []
    startmin = 0
    while True:
        endmin = min(startmin + 1, duration)
        intervals.append((startmin, endmin))
        if endmin == duration:
            break
        startmin = endmin

    t0 = time.time()
    l = []
    for start, end in intervals:
        raw = load_ephys(path, startmin=start, endmin=end, offset=lfpstart, cleaned=True)
        X = raw_lfp_to_csd(raw.values)
        idx = raw.index[::2]
        col = raw.columns[::4] // 4
        csd = pd.DataFrame(X, index=idx, columns=col)
        l.append(csd)

    csd = pd.concat(l, axis=0)
    logging.info(f'total samples: {len(csd)}')
    logging.info(f'GB: {csd.values.nbytes / 1e9}')
    return csd


def upsample_behavior(date):
    with open(src.paths.REF / 'date_recid.json', 'r') as f:
        recids = json.load(f)
    with open(src.paths.REF / 'date_recoffset.json', 'r') as f:
        offsets = json.load(f)
        
    recid = recids[date]
    lfpstart, lfpend = offsets[date]

    df = load_behavior(recid)
    lfplen = (lfpend - lfpstart + 1) // 2
    lfp_t = pd.date_range(start=0, periods=lfplen, freq='0.4ms')
    behav_t = pd.date_range(start=0, periods=len(df), freq='10ms')

    df = df.set_index(behav_t)
    df = df.reindex(lfp_t).interpolate(method='linear')
    df.index = pd.Series(np.arange(len(df)) / 2500, name='time')
    return df
