## Setup

import os, sys
from pathlib import Path
sys.path.append('/Users/Admin/Documents/scripts/effie')
sys.path.append('/Users/Admin/Documents/scripts/trigrid')


import cartopy.crs as ccrs, pandas as pd

from effie.grid import Grid
import effie.source, effie.target, effie.interface

top_level_str = '/Users/Admin/Documents/scripts/fvcom-work'

## Grid

g = Grid.from_fvcom_dat(grd_file=   os.path.join(top_level_str, 'Adamselv/input/adamselv_grd.dat'),
                        obc_file=   os.path.join(top_level_str, 'Adamselv/input/adamselv_obc.dat'),
                        dep_file=   os.path.join(top_level_str, 'Adamselv/input/adamselv_dep.dat'),
                        sigma_descr=os.path.join(top_level_str, 'Adamselv/input/adamselv_sigma.dat'),
                        crs=ccrs.UTM(33))


## Input files

t0 = pd.Timestamp('2018-2-1')
t1 = pd.Timestamp('2018-6-30')
outdir = os.path.join(top_level_str, 'Adamselv/03_atm_forcing/output')


## Atmospheric forcing

source = effie.source.meps.MEPS(start_time=t0, end_time=t1, cache=True)
target = effie.target.forcing.Forcing(grid=g, casename='adamselv', outdir=outdir)
i = effie.interface.NdInterpolation(source=source, target=target)
i.write_target_file(clobber=True)
