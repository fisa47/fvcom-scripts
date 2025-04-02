import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegWriter
import contextily as ctx
import pandas as pd
import numpy as np

# Load dataset
ds_all = xr.open_dataset('../run-output/pipe_ave/adamselv_v01_0001.nc', decode_times=False)


def correct_tracer_concentration(ds, tracer_name, vol, time, pipe_concentration=10, pipe_discharge=5.3/60):
    """
    Applies mass correction to a tracer and returns corrected concentration.

    Parameters:
    -----------
    ds : xarray.Dataset
        Dataset containing the tracer.
    tracer_name : str
        Name of the tracer variable in the dataset (e.g. 'tracer1_c').
    vol : np.ndarray
        3D array of volume per cell [time, layers, elements].
    time : np.ndarray
        1D array of days since FVCOM origin
    pipe_concentration : float
        Tracer concentration at the pipe (default 10).
    pipe_discharge : float
        Discharge rate from the pipe in m³/s (default 5.3/60).
    
    Returns:
    --------
    tracer_conc : np.ndarray
        Corrected concentration of the tracer [time, layers, elements].
    """
    tracer_mass = vol * ds[tracer_name].values
    tracer_mass_discharged = pipe_concentration * pipe_discharge * (np.abs(time[0] - time) * 86400)

    int_mass = tracer_mass.sum(axis=(-1, -2))
    calibration = tracer_mass_discharged / int_mass
    calibration[0] = 1  # avoid NaN or inf in the first timestep

    corrected_mass = tracer_mass * calibration[:, np.newaxis, np.newaxis]
    tracer_conc = corrected_mass / vol
    return tracer_conc


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

print(ds_all)

tracer1_conc = correct_tracer_concentration(ds_all, 'tracer1_c', vol, time)
tracer2_conc = correct_tracer_concentration(ds_all, 'tracer2_c', vol, time)

# Get the 2nd from the bottom layer
ds_all = ds_all.isel(siglay=1)
# Apply correction to both tracers and NB!!! get max concentration across depth layers
ds_all['tracer1_c'].values = tracer1_conc.max(axis=1)  # max across depth layers
ds_all['tracer2_c'].values = tracer2_conc.max(axis=1)


# Global min/max for consistent color scaling
vminmax = {
    var: (ds_all[var].min().item(), ds_all[var].max().item())
    for var in ['temp', 'salinity', 'tracer1_c', 'tracer2_c', 'ua', 'va']
}

vminmax['tracer1_c'] = (0, 1)

# Setup figure and axes
fig, axs = plt.subplots(3, 2, figsize=(12, 12), sharex=True, sharey=True)
plots = {}
colorbars = {}

time_discharge = time[28] # 14th day

def init():
    ds = ds_all.isel(time=0)
    variables = ['temp', 'salinity', 'tracer1_c', 'tracer2_c', 'ua', 'va']
    titles = ['temperature', 'salinity', 'tracer1_c, MAX in water column', 'tracer2_c, MAX in water column', 'u', 'v']
    cmaps = ['plasma', 'viridis', 'magma_r', 'magma_r', 'PuBuGn', 'PuBuGn']
    positions = [(0, 0), (0, 1), (1, 0), (1, 1), (2, 0), (2, 1)]

    for var, title, cmap, pos in zip(variables, titles, cmaps, positions):
        ax = axs[pos]
        if var in ['ua', 'va']:
            x, y = ds['xc'], ds['yc']
        else:
            x, y = ds['x'], ds['y']

        plots[var] = ax.scatter(x, y, c=ds[var], cmap=cmap, s=1,
                                vmin=vminmax[var][0], vmax=vminmax[var][1])
        ax.set_title(title)
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
    current_time = time[i]
    # time is in days since FVCOM origin
    delta_days = (current_time - time_discharge)

    for var in plots:
        plots[var].set_array(ds[var].values)

    if delta_days >= 0:
        fig.suptitle(f"Days since discharge: {delta_days:.2f}", fontsize=16)
    else:
        fig.suptitle(f"Days before discharge: {abs(delta_days):.2f}", fontsize=16)

    print(f"Frame {i} updated")
    return plots.values()

# Animate
anim = FuncAnimation(fig, update, init_func=init, frames=60, blit=False)

# Save as MP4
writer = FFMpegWriter(fps=1, metadata=dict(artist='Python'), bitrate=1800)
anim.save("output_animation3.mp4", writer=writer)

plt.close()
