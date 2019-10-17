#!/usr/env/python

## Import General Tools
import sys
import os
import argparse
import logging
from datetime import datetime as dt
from datetime import timedelta as tdelta

from pathlib import Path
import subprocess

import numpy as np
from astropy import units as u
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS, FITSFixedWarning
from astropy.table import QTable, Row

import warnings
warnings.simplefilter('ignore', FITSFixedWarning)

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
p.add_argument("--relax", dest="relax",
    default=False, action="store_true",
    help="Don't verify pointing origin and other keywords (for testing).")
## add options
# p.add_argument("--input", dest="input", type=str,
#     help="The input.")
## add arguments
p.add_argument('files', nargs='*',
               help="All files")
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
LogFileName = Path('ImageResults.log').expanduser()
LogFileHandler = logging.FileHandler(LogFileName)
LogFileHandler.setLevel(logging.DEBUG)
LogFileHandler.setFormatter(LogFormat)
log.addHandler(LogFileHandler)

##-------------------------------------------------------------------------
## Get Center Coord
##-------------------------------------------------------------------------
def get_center_coord(hdu):
    w = WCS(hdu.header)
    nx, ny = hdu.data.shape
    result = w.all_pix2world(np.array(nx/2), np.array(ny/2), 1)
    center_coord = SkyCoord(result[0], result[1], frame='icrs', unit='deg')
    return(center_coord)

##-------------------------------------------------------------------------
## solve_pointing
##-------------------------------------------------------------------------
def solve_pointing(filename, relax=False):
    tick = dt.utcnow()
    log.info(f"Analyzing {filename.name}")
    hdr = fits.getheader(filename)
    if not relax:
        if not hdr.get('PONAME').strip() == 'REF': return None
        if not float(hdr.get('RAOFF')) < 0.1: return None
        if not float(hdr.get('DECOFF')) < 0.1: return None
        if not hdr.get('OBJECT') in ['GuiderFlexureTest', 'MIRA PMFM 350', 'MIRA PMFM -350']: return None

    # Extract EL, SKYPA, ROTPPOSN, FILTER
    EL = float(hdr.get('EL')) * u.deg
    SKYPA1 = float(hdr.get('SKYPA1')) * u.deg
    SKYPA2 = float(hdr.get('SKYPA2')) * u.deg
    ROTPPOSN = float(hdr.get('ROTPPOSN')) * u.deg
    FILTER = hdr.get('FILTER')
    OBJECT = hdr.get('OBJECT')
    FCPA, FCEL = (hdr.get('FCPA_EL')).split(' ')

    # Solve image for astrometry
    astrometry_cmd = ['solve-field', '-O', '-p', '-z', '2', '-t', '2', f"{f}"]
#     astrometry_cmd = ['solve-field', '-O', '-p', '-z', '2', '-T', f"{f}"]
    log.info('  ' + ' '.join(astrometry_cmd[:-1]))
    output = subprocess.run(astrometry_cmd, stdout=subprocess.PIPE)
    solved_file = f.with_name(f.name.replace('.fits', '.solved'))
    new_file = f.with_name(f.name.replace('.fits', '.new'))
    if not solved_file.exists():
        log.error(f'  Astrometry solve failed for {filename.name}')
        return
    # Open Solved Image
    hdul = fits.open(new_file)

    # Extract Guider Pointing Info from Header
    PONAME = hdul[0].header.get('PONAME') # This should be REF
    POYPOS = hdul[0].header.get('POYPOS')
    POXPOS = hdul[0].header.get('POXPOS')
    POYOFF = hdul[0].header.get('POYOFF')
    POXOFF = hdul[0].header.get('POXOFF')
    TARGRA = hdul[0].header.get('TARGRA')
    TARGDE = hdul[0].header.get('TARGDE')
    TARGFR = hdul[0].header.get('TARGFR')
    ROTSTS = hdul[0].header.get('ROTSTS')
    ROTMOD = hdul[0].header.get('ROTMOD')
    RAOFF = float(hdul[0].header.get('RAOFF')) * u.arcsec
    DECOFF = float(hdul[0].header.get('DECOFF')) * u.arcsec
    RA = float(hdul[0].header.get('RA')) * u.deg
    DEC = float(hdul[0].header.get('DEC')) * u.deg

    guider_coord = SkyCoord(RA, DEC, frame='icrs')

    # Extract WCS from header
    center_coord = get_center_coord(hdul[0])
    
#     w = WCS(hdul[0].header)
#     nx, ny = hdul[0].data.shape
#     result = w.all_pix2world(np.array(nx/2), np.array(ny/2), 1)
#     center_coord = SkyCoord(result[0], result[1], frame='icrs', unit='deg')

    offset_pa = center_coord.position_angle(guider_coord).to(u.deg)
    offset_distance = center_coord.separation(guider_coord).to(u.arcsec)
    offset_angle = offset_pa.to(u.deg) - SKYPA2
    offset_angle.wrap_at(180*u.deg, inplace=True)

    tock = dt.utcnow()
    analysis_time = (tock-tick).total_seconds()
    log.info(f"  Solved {analysis_time:.0f} s: {offset_distance:.1f}, "\
             f"{offset_angle:.2f} at drive = {ROTPPOSN:.2f}")

    table_file = Path('ImageResults.txt').expanduser()
    if not table_file.exists():
        t = QTable()
        t['Filename'] = [f.name]
        t['EL'] = [EL.value]
        t['PA'] = [SKYPA2.value]
        t['RotAng'] = [ROTPPOSN.value]
        t['GuiderCoord'] = [guider_coord.to_string(style='hmsdms', sep=':', precision=2)]
        t['ImageCoord'] = [center_coord.to_string(style='hmsdms', sep=':', precision=2)]
        t['OffsetDistance'] = [offset_distance.value]
        t['OffsetAngle'] = [offset_angle.value]
        print(t)
    else:
        t = QTable.read(table_file, format='ascii.ecsv')
        row = {'Filename': f.name,
               'EL': EL.value,
               'PA': SKYPA2.value,
               'RotAng': ROTPPOSN.value,
               'GuiderCoord': guider_coord.to_string(style='hmsdms', sep=':', precision=2),
               'ImageCoord': center_coord.to_string(style='hmsdms', sep=':', precision=2),
               'OffsetDistance': offset_distance.value,
               'OffsetAngle': offset_angle.value,
              }
        t.add_row(vals=row)
    t.write(table_file, format='ascii.ecsv', overwrite=True)


if __name__ == '__main__':
   
    for file in args.files:
        f = Path(file).expanduser()
        solve_pointing(f, relax=args.relax)
