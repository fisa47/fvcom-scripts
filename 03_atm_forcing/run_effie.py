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

t0 = pd.Timestamp('2017-2-1')
t1 = pd.Timestamp('2017-6-30')
outdir = os.path.join(top_level_str, 'Adamselv/atm_forcing/output')


## Atmospheric forcing

source = effie.source.meps.MEPS(start_time=t0, end_time=t1, cache=True)
target = effie.target.forcing.Forcing(grid=g, casename='adamselv', outdir=outdir)
i = effie.interface.NdInterpolation(source=source, target=target)
i.write_target_file(clobber=True)


# ### Boundary (nesting) file

# source = effie.source.roms.NorKyst800(start_time=t0, end_time=t1, gridpoints_per_tile=100**2, cache=True)
# target = effie.target.boundary.BoundaryType3Nesting(
#     grid=g, casename='case', nesting_zone_width_layers=4, outdir=outdir,
#     ramp={v: (val, t0+pd.Timedelta(hours=1), t0+pd.Timedelta(days=2)) for v, val in zip(
#         ['u', 'v', 'zeta', 'temp', 'salinity'],
#         [0, 0, 0, 6, 35],
#     )}
# )
# i = effie.interface.NdInterpolation(source=source, target=target)
# i.write_target_file(clobber=True)


# ### Initial (restart) file

# source = effie.source.roms.NorKyst800Snapshot(time=t0, gridpoints_per_tile=100**2, cache=True)
# target = effie.target.state.Initial(grid=g, casename='case', outdir=outdir)
# i = effie.interface.NdInterpolation(source=source, target=target)
# i.write_target_file(clobber=True)
