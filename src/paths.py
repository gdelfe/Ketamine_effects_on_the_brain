import glob
import os
import pathlib


ROOT = pathlib.Path(__file__).parent.parent
DATA = ROOT / 'data'
REF = ROOT / 'ref'


def session_path(date):
    datadir = DATA / date
    if not os.path.isdir(datadir):
        raise ValueError(f"data directory {datadir} not found")
    
    # get the path from a file that's already there
    pattern = datadir / f'*.imec1.lf.meta'
    match = glob.glob(str(datadir / pattern))[0]
    file = os.path.split(match)[-1]
    return file.replace(f'.imec1.lf.meta', '')
