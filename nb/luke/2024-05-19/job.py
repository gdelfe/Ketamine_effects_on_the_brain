import argparse
import json
import logging
import os
import pickle
import sys

from src.pipeline import compute_session_csd


def main(output_path, params):
    pattern = '_'.join([f'{k}={v}' for k, v in params.items()])
    outfile = os.path.join(output_path, pattern + '.pkl')
    if os.path.isfile(outfile):
        logging.info(f'output file {outfile} already exists, exiting')
        sys.exit(0)
    
    logging.info(f'computing csd for recording {params}')
    result = compute_session_csd(params['date'], params['region'])
    logging.info('done computing csd')
    
    logging.info(f'writing result to {outfile}')
    with open(outfile, 'wb') as f:
        pickle.dump(result, f)
    logging.info('done')

    
if __name__ == '__main__':
    logging.basicConfig(
        level=logging.INFO,
        format=f'[%(asctime)s] %(message)s',
        handlers=[logging.StreamHandler()]
    )
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--output-path')
    parser.add_argument('-p', '--params')
    args = parser.parse_args()
    main(args.output_path, json.loads(args.params))
