#!/usr/env/python

## Import General Tools
import sys
import os
import argparse
import logging

from pathlib import Path
import subprocess

import numpy as np
from astropy import units as u
from astropy.io import fits
from astropy import coordinates as c
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
## solve_pointing
##-------------------------------------------------------------------------
def solve_pointing(filename, relax=False):
    

    # Solve image for astrometry
    output = subprocess.run(['solve-field', '-O', '-p', f"{f}"],
                            stdout=subprocess.PIPE)
    solved_file = f.with_name(f.name.replace('.fits', '.solved'))
    new_file = f.with_name(f.name.replace('.fits', '.new'))
    if not solved_file.exists():
        print('Astrometry solve failed')
        raise FileNotFoundError

    # Open Solved Image
    hdul = fits.open(new_file)

    # Extract EL, SKYPA, ROTPPOSN, FILTER
    EL = float(hdul[0].header.get('EL')) * u.deg
    SKYPA2 = float(hdul[0].header.get('SKYPA2')) * u.deg
    ROTPPOSN = float(hdul[0].header.get('ROTPPOSN')) * u.deg
    FILTER = hdul[0].header.get('FILTER')

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
    if not relax:
        assert PONAME.strip() == 'REF'
        assert RAOFF < 0.1 * u.arcsec
        assert DECOFF < 0.1 * u.arcsec

    guider_coord = c.SkyCoord(RA, DEC, frame='icrs')

    # Extract WCS from header
    w = WCS(hdul[0].header)
    nx, ny = hdul[0].data.shape
    result = w.all_pix2world(np.array(nx/2), np.array(ny/2), 0)
    center_coord = c.SkyCoord(result[0], result[1], frame='icrs', unit='deg')

    offset_pa = center_coord.position_angle(guider_coord).to(u.deg)
    offset_distance = center_coord.separation(guider_coord).to(u.arcsec)
    physical_offset_angle = offset_pa.to(u.deg) - SKYPA2

    print(f"{f.name:18s}: {offset_distance:.1f} at PA = {offset_pa:.2f}, SKYPA2={SKYPA2:.2f} (difference={physical_offset_angle:.2f})")

    table_file = Path('~/KeckData/MOSFIRE_GuiderFlexure/ImageResults.txt').expanduser()
    if not table_file.exists():
        t = QTable()
        t['Filename'] = [f.name]
        t['EL'] = [EL.value] * EL.unit
        t['PA'] = [SKYPA2.value] * SKYPA2.unit
        t['RotAng'] = [ROTPPOSN.value] * ROTPPOSN.unit
        t['GuiderCoord'] = [guider_coord.to_string(style='hmsdms', sep=':', precision=2)]
        t['ImageCoord'] = [center_coord.to_string(style='hmsdms', sep=':', precision=2)]
        t['Offset Distance'] = [offset_distance.value] * offset_distance.unit
        t['Offset Angle'] = [physical_offset_angle.value] * physical_offset_angle.unit
        print(t)
    else:
        t = QTable.read(table_file, format='ascii')
        row = {'Filename': f.name,
               'EL': EL,
               'PA': SKYPA2,
               'RotAng': ROTPPOSN,
               'GuiderCoord': guider_coord.to_string(style='hmsdms', sep=':', precision=2),
               'ImageCoord': center_coord.to_string(style='hmsdms', sep=':', precision=2),
               'Offset Distance': offset_distance,
               'Offset Angle': physical_offset_angle,
              }
        print(row)
        t.add_row(vals=row)

    t.write(table_file, format='ascii.ecsv', overwrite=True)

if __name__ == '__main__':
   
    for file in args.files:
        f = Path(file).expanduser()
        solve_pointing(f, relax=args.relax)
