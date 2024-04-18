import os
import pandas as pd

sessions = pd.read_csv('/Users/lukearend/phd/kfx/ref/sessions.csv', index_col=0)
recordings = pd.read_csv('/Users/lukearend/phd/kfx/ref/recordings.csv', index_col=0)
metadata = pd.read_csv('/Users/lukearend/phd/kfx/ref/metadata.csv', index_col=[0, 1, 2])


def data_path(sesspath, signal='lfp', region='hpc', ftype='data'):
    probe = {'pfc': 'imec0', 'hpc': 'imec1'}[region.lower()]
    sig = {'lfp': 'lf', 'ap': 'ap'}[signal.lower()]
    ext = {'data': 'bin', 'meta': 'meta'}[ftype]
    sessdir, sesspath = os.path.split(sesspath)
    return f'{sessdir}/{sesspath}_{probe}/{sesspath}_t0.{probe}.{sig}.{ext}'


def load_metadata(path):
    with open(path, 'r') as f:
        metadata = {}
        for line in f.read().splitlines():
            k, v = line.split('=')
            k = k.strip('~')
            metadata[k] = v
    return metadata