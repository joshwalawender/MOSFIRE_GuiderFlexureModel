#!/usr/env/python

## Import General Tools
import sys
import os
import argparse
import logging

from pathlib import Path
import subprocess

import numpy as np
from astropy.io import fits
from astropy import coordinates as c
from astropy.wcs import WCS

##-------------------------------------------------------------------------
## Parse Command Line Arguments
##-------------------------------------------------------------------------
## create a parser object for understanding command-line arguments
p = argparse.ArgumentParser(description='''
''')
## add flags
p.add_argument("-v", "--verbose", dest="verbose",
    default=False, action="store_true",
    help="Be verbose! (default = False)")
## add options
# p.add_argument("--input", dest="input", type=str,
#     help="The input.")
## add arguments
p.add_argument('file', type=str,
               help="A single argument")
p.add_argument('allothers', nargs='*',
               help="All other arguments")
args = p.parse_args()


##-------------------------------------------------------------------------
## Create logger object
##-------------------------------------------------------------------------
log = logging.getLogger('MyLogger')
log.setLevel(logging.DEBUG)
## Set up console output
LogConsoleHandler = logging.StreamHandler()
LogConsoleHandler.setLevel(logging.DEBUG)
LogFormat = logging.Formatter('%(asctime)s %(levelname)8s: %(message)s',
                              datefmt='%Y-%m-%d %H:%M:%S')
LogConsoleHandler.setFormatter(LogFormat)
log.addHandler(LogConsoleHandler)
## Set up file output
# LogFileName = None
# LogFileHandler = logging.FileHandler(LogFileName)
# LogFileHandler.setLevel(logging.DEBUG)
# LogFileHandler.setFormatter(LogFormat)
# log.addHandler(LogFileHandler)



##-------------------------------------------------------------------------
## solve_pointing
##-------------------------------------------------------------------------
def solve_pointing(filename):
    f = Path(filename).expanduser()
    
    if not f.exists(): raise FileNotFoundError

    subprocess.call(['solve-field', '-O', '-p', f"{f}"])
    solved_file = f.with_name(f.name.replace('.fits', '.solved'))
    new_file = f.with_name(f.name.replace('.fits', '.new'))
    if not solved_file.exists():
        print('Astrometry solve failed')
        raise FileNotFoundError


    hdul = fits.open(new_file)
    el = hdul[0].header.get('EL')
    skypa2 = hdul[0].header.get('SKYPA2')
    rotpposn = hdul[0].header.get('ROTPPOSN')
    filt = hdul[0].header.get('FILTER')

    w = WCS(hdul[0].header)
    nx, ny = hdul[0].data.shape
    result = w.all_pix2world(np.array(nx/2), np.array(ny/2), 0)
    center_coord = c.SkyCoord(result[0], result[1], frame='icrs', unit='deg')

    print(f"{f.name:20s}: {el:.1f} {skypa2:+6.1f} {rotpposn:+6.1f} "\
          f"{center_coord.ra.deg:.5f} {center_coord.dec.deg:.5f}")

if __name__ == '__main__':
    solve_pointing(args.file)
