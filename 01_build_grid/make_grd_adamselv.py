import numpy as np
import datetime as dt
import netCDF4 as nc
import sys
sys.path.append('/Users/Admin/Documents/scripts/pyfvcom')

import matplotlib
font = {'size'   : 16}

matplotlib.rc('font', **font)

from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.pyplot as plt

import PyFVCOM as pf

from pyproj import Proj
utm35 = Proj(proj = 'utm', zone = 35, ellps = 'WGS84', preserve_units = False)


def add_cbar(ax, fig, tc, label):
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    cbar = fig.colorbar(tc, cax=cax, orientation='vertical')
    cbar.set_label(label)


project = 'Adamselv'
casename = 'adamselv_v01'
grid_2dm = f'{casename}.2dm'

parent_case = 'finntest' 
parent_mesh = f'input/{parent_case}_grd.dat'
parent_dep = f'input/{parent_case}_dep.dat'

tri, nodes, x, y, z, types, obc = pf.grid.read_sms_mesh(grid_2dm, nodestrings=True)

old_tri, _, old_x, old_y, _ = pf.grid.read_fvcom_mesh(parent_mesh)
old_dep = np.loadtxt(parent_dep, skiprows=1)


choose_poly = np.array([[26.6673 , 70.9506],
       [27.2379 , 70.9447 ],
       [27.8004 , 70.8797 ],
       [27.1447 , 70.2766 ],
       [26.4987 , 70.2938 ],
       [26.20853, 70.60794],
       [26.6673 , 70.9506],])


import shapely.geometry as sg
subset_poly = sg.Polygon(choose_poly)
subset_poly = subset_poly.buffer(5000)

choose_dep = [subset_poly.contains(sg.Point(x,y)) for x,y in zip(old_dep[:,0], old_dep[:,1])]
print('sum', np.sum(choose_dep))  # Check how many points are selected


from scipy.interpolate import LinearNDInterpolator, NearestNDInterpolator, griddata

lon, lat = utm35(x,y, inverse=True)

ldi = LinearNDInterpolator(old_dep[choose_dep, 0:2], old_dep[choose_dep,2])
z_old = ldi(np.vstack([lon,lat]).T)
nni = NearestNDInterpolator(old_dep[choose_dep, 0:2], old_dep[choose_dep,2])
z_old = nni(np.vstack([lon,lat]).T)

print(z_old.max(), z_old.min())

print('xy', old_dep[choose_dep, :2])  # Should contain valid (x, y) pairs
print('depth', old_dep[choose_dep, 2])   # Should contain depth values

print('nans', np.isnan(z_old).sum())  # How many NaNs exist before nearest-neighbor interpolation?


# EMODnet 

lon, lat = utm35(x,y, inverse=True)
emod_nc = nc.Dataset('C7_2022.nc')

offset = 0.2
choose_lon = np.logical_and(emod_nc['lon'][:] >= np.min(lon) - offset, emod_nc['lon'][:] <= np.max(lon) + offset)
choose_lat = np.logical_and(emod_nc['lat'][:] >= np.min(lat) - offset, emod_nc['lat'][:] <= np.max(lat) + offset)

plt.pcolormesh(emod_nc['lon'][choose_lon], emod_nc['lat'][choose_lat], emod_nc['elevation'][choose_lat, :][:, choose_lon])
plt.triplot(lon,lat,tri,c='r', alpha=0.5, zorder=2)

emod_ll = np.meshgrid(emod_nc['lon'][choose_lon], emod_nc['lat'][choose_lat])
emod_xy = utm35(emod_ll[0], emod_ll[1])
emod_elev = emod_nc['elevation'][choose_lat, :][:, choose_lon]



ldi = LinearNDInterpolator(np.asarray([emod_xy[0].ravel(), emod_xy[1].ravel()]).T, emod_elev.ravel())
z_emod = ldi(np.vstack([x,y]).T)

rem = np.isnan(emod_elev.ravel())

nni = NearestNDInterpolator(np.asarray([emod_xy[0].ravel()[~rem], emod_xy[1].ravel()[~rem]]).T, emod_elev.ravel()[~rem])
z_emod[np.isnan(z_emod)] = nni(np.vstack([x,y]).T[np.isnan(z_emod),:])

to_plot = {
    'zoom_out': {
        'x_limits': [470804, 511239,],
        'y_limits': [7807577, 7884528],
        'depth_range': [0, 300]
    },
    'zoom_in': {
        'x_limits': [486067, 493642],
        'y_limits': [7811617, 7823099],
        'depth_range': [0, 200]
    }
}

for name, limits in to_plot.items():
    fig, axs = plt.subplots(1, 2, figsize=(18, 10))
    
    # Plot old bathymetry
    tc = axs[0].tripcolor(x, y, tri, z_old, 
                          vmax=limits['depth_range'][1], 
                          vmin=limits['depth_range'][0])
    axs[0].set_xlim(limits['x_limits'])
    axs[0].set_ylim(limits['y_limits'])
    axs[0].set_aspect('equal')
    add_cbar(axs[0], fig, tc, 'Depth (m)')

    # Plot EMODnet bathymetry
    tc = axs[1].tripcolor(x, y, tri, -z_emod, 
                          vmax=limits['depth_range'][1], 
                          vmin=limits['depth_range'][0])
    axs[1].set_xlim(limits['x_limits'])
    axs[1].set_ylim(limits['y_limits'])
    axs[1].set_aspect('equal')
    add_cbar(axs[1], fig, tc, 'Depth (m)')

    fig.tight_layout()
    fig.savefig(f'bathy_{name}.png', dpi=180)
    plt.close()



"""
# dybdata - can't find where to download this anymore...

xyz = np.loadtxt('xyz_0.xyz')

ldi = LinearNDInterpolator(xyz[:,0:2], xyz[:,2])
z = ldi(np.vstack([x,y]).T)
nni = NearestNDInterpolator(xyz[:,0:2], xyz[:,2])
z[np.isnan(z)] = nni(np.vstack([x,y]).T[np.isnan(z),:])
"""

# patch together emodnet for round the islands with the old bathymetry outside, smoothing should be done by fvtools preproc routines assuming the match isn't really terrible

choose_dyb = np.array([[26.66 , 70.406],
       [26.804 , 70.412 ],
       [26.75, 70.488],
       [26.626 , 70.488],
       [26.66 , 70.406],
       ])

dyb_use_poly = sg.Polygon(choose_dyb)
choose_dyb = [dyb_use_poly.contains(sg.Point(x1,y1)) for x1,y1 in zip(x,y)]

hybrid_z = z_old
hybrid_z[choose_dyb] = -z_emod[choose_dyb]

np.save(f'bathy_hybrid_{casename}.npy', np.asarray([x,y,hybrid_z]).T)

"""
pf.grid.write_fvcom_mesh(tri, nodes, x, y, z, f'{casename}_grd.dat', extra_depth=f'{casename}_dep.dat')
pf.grid.write_obc_file(np.squeeze(obc), np.ones(len(np.squeeze(obc))), f'{casename}_obc.dat')

model = pf.preproc.Model(start=dt.datetime(2020,1,1), end=dt.datetime(2020,4,1), grid=f'{casename}_grd.dat', native_coordinates='cartesian', zone='33')
model.write_coriolis(f'{casename}_cor.dat')
"""



