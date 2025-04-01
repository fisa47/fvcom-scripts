import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegWriter
import contextily as ctx
import pandas as pd
import numpy as np

# Load dataset
ds_all = xr.open_dataset('../run-output/pipe_ave/adamselv_v01_0001.nc', decode_times=False)


# Volume preprocessing
def unstructured_grid_volume(area, depth, surface_elevation, thickness):
    dz = np.abs(np.diff(thickness, axis=0))
    volume = (area * (surface_elevation + depth))
    depth_volume = volume[:, np.newaxis, :] * dz[np.newaxis, ...]
    return depth_volume

# Load required vars
area = ds_all['art1'].values
depth = ds_all['h'].values
surface_elevation = ds_all['zeta'].values
thickness = ds_all['siglev'].values
time = ds_all.time.values

# Calculate volume and mass
vol = unstructured_grid_volume(area, depth, surface_elevation, thickness)
tracer1_mass = vol * ds_all['tracer1_c'].values
tracer2_mass = vol * ds_all['tracer2_c'].values

# Pipe discharge setup (adjust as needed)
pipe_c = 10
pipe_discharge = 5.3 / 60  # m3/s
tracer_mass_discharged = pipe_c * pipe_discharge * (np.abs(time[0] - time) * 86400) 

# Mass conservation correction
int_mass1 = tracer1_mass.sum(axis=(-1, -2))
int_mass2 = tracer2_mass.sum(axis=(-1, -2))

cal1 = tracer_mass_discharged / int_mass1
cal2 = tracer_mass_discharged / int_mass2
cal1[0] = 1
cal2[0] = 1

# Apply correction
tracer1_mass = tracer1_mass * cal1[:, np.newaxis, np.newaxis]
tracer2_mass = tracer2_mass * cal2[:, np.newaxis, np.newaxis]

# Recalculate concentration
tracer1_conc = tracer1_mass / vol
tracer2_conc = tracer2_mass / vol

# Assign to dataset
ds_all['tracer1_c'].values = tracer1_conc
ds_all['tracer2_c'].values = tracer2_conc

print(ds_all)

ds_all = ds_all.isel(siglay=1)

# Global min/max for consistent color scaling
vminmax = {
    var: (ds_all[var].min().item(), ds_all[var].max().item())
    for var in ['temp', 'salinity', 'tracer1_c', 'tracer2_c', 'ua', 'va']
}

# Setup figure and axes
fig, axs = plt.subplots(3, 2, figsize=(15, 15), sharex=True, sharey=True)
plots = {}
colorbars = {}

def init():
    ds = ds_all.isel(time=0)
    variables = ['temp', 'salinity', 'tracer1_c', 'tracer2_c', 'ua', 'va']
    cmaps = ['plasma', 'viridis', 'magma_r', 'magma_r', 'PuBuGn', 'PuBuGn']
    positions = [(0, 0), (0, 1), (1, 0), (1, 1), (2, 0), (2, 1)]

    for var, cmap, pos in zip(variables, cmaps, positions):
        ax = axs[pos]
        if var in ['ua', 'va']:
            x, y = ds['xc'], ds['yc']
        else:
            x, y = ds['x'], ds['y']

        plots[var] = ax.scatter(x, y, c=ds[var], cmap=cmap, s=1,
                                vmin=vminmax[var][0], vmax=vminmax[var][1])
        ax.set_title(var)
        ax.set_xlim(9.32e5, 9.38e5)
        ax.set_ylim(7.853e6, 7.86e6)
        ax.grid(True)

        # Add basemap
        ctx.add_basemap(ax, crs='EPSG:32633', source=ctx.providers.CartoDB.Positron, attribution_size=6)

        # Add static colorbar
        colorbars[var] = fig.colorbar(plots[var], ax=ax, label=var)

    return plots.values()

def update(i):
    ds = ds_all.isel(time=i)
    plots['temp'].set_array(ds['temp'].values)
    plots['salinity'].set_array(ds['salinity'].values)
    plots['tracer1_c'].set_array(ds['tracer1_c'].values)
    plots['tracer2_c'].set_array(ds['tracer2_c'].values)
    plots['ua'].set_array(ds['ua'].values)
    plots['va'].set_array(ds['va'].values)
    print(f"Frame {i} updated")
    return plots.values()

# Animate
anim = FuncAnimation(fig, update, init_func=init, frames=60, blit=False)

# Save as MP4
writer = FFMpegWriter(fps=1, metadata=dict(artist='Python'), bitrate=1800)
anim.save("output_animation2.mp4", writer=writer)

plt.close()
