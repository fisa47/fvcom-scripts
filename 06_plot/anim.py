import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegWriter
import contextily as ctx
import pandas as pd
import numpy as np
from datetime import datetime, timedelta

# FVCOM origin date
fvcom_origin = datetime(1858, 11, 17)

pipe_discharge = 5.3/60  # m3/s, from Grieg email

# Nutrients data from
# Notat 62802, Sak: Søknad om endring i utslippstillatelse Adamselv Settefisk
# Table 5
yr_vol = pipe_discharge * 86400 * 365  # m3
N_g_m3 = 186305*1000/yr_vol  # 66.9 g/m3
P_g_m3 = 28622*1000/yr_vol  # 10.3  g/m3
TOC_g_m3 = 342054*1000/yr_vol  # 122.8 g/m3
pipe_ind = 54469

# Load dataset
ds_all = xr.open_mfdataset('../run-output/pipe_ave/adamselv_v01_000*.nc', decode_times=False, concat_dim='time', combine='nested')

def correct_tracer_concentration(ds, tracer_name, vol, time,
                                 pipe_concentration=10, pipe_discharge=5.3/60,
                                 substance_concentration=10):
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
    # Step 2. Correct for the mass of the tracer given pipe_concentration in the model (10 g/m3)
    tracer_mass = vol * ds[tracer_name].values
    tracer_mass_discharged = pipe_concentration * pipe_discharge * (np.abs(time[0] - time) * 86400)

    # model mass released at each timestep
    int_mass = tracer_mass.sum(axis=(-1, -2))
    calibration = tracer_mass_discharged / int_mass
    calibration[0] = 1  # avoid NaN or inf in the first timestep

    corrected_mass = tracer_mass * calibration[:, np.newaxis, np.newaxis]

    # Step 2. Recalculate for the mass of the substance
    substance_mass_discharged = substance_concentration * pipe_discharge * (np.abs(time[0] - time) * 86400)

    # mass released at each timestep
    int_corrected_mass = corrected_mass.sum(axis=(-1, -2))
    calibration_substance = substance_mass_discharged  / int_corrected_mass
    corrected_substance_mass = corrected_mass * calibration_substance[:, np.newaxis, np.newaxis]

    tracer_conc = corrected_substance_mass / vol
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
tri = ds_all['nv'].isel(time=0).values.T - 1


# Calculate volume and mass
vol = unstructured_grid_volume(area, depth, surface_elevation, thickness)

# tracer1_conc = correct_tracer_concentration(ds_all, 'tracer1_c', vol, time)
# tracer2_conc = correct_tracer_concentration(ds_all, 'tracer2_c', vol, time)

N_conc = correct_tracer_concentration(ds_all, 'tracer1_c', vol, time, substance_concentration=N_g_m3)
N_conc_heavy = correct_tracer_concentration(ds_all, 'tracer2_c', vol, time, substance_concentration=N_g_m3)

P_conc = correct_tracer_concentration(ds_all, 'tracer1_c', vol, time, substance_concentration=P_g_m3)
P_conc_heavy = correct_tracer_concentration(ds_all, 'tracer2_c', vol, time, substance_concentration=P_g_m3)

TOC_conc = correct_tracer_concentration(ds_all, 'tracer1_c', vol, time, substance_concentration=TOC_g_m3)
TOC_conc_heavy = correct_tracer_concentration(ds_all, 'tracer2_c', vol, time, substance_concentration=TOC_g_m3)


# Create new variables for corrected concentrations
ds_all['N_c'] = ds_all['tracer1_c'] * 0
ds_all['N_c_heavy'] = ds_all['tracer1_c'] * 0

ds_all['P_c'] = ds_all['tracer1_c'] * 0
ds_all['P_c_heavy'] = ds_all['tracer1_c'] * 0

ds_all['TOC_c'] = ds_all['tracer1_c'] * 0
ds_all['TOC_c_heavy'] = ds_all['tracer1_c'] * 0

# Assign corrected concentrations to new variables
ds_all['N_c'].values = N_conc
ds_all['N_c_heavy'].values = N_conc_heavy

ds_all['P_c'].values = P_conc
ds_all['P_c_heavy'].values = P_conc_heavy

ds_all['TOC_c'].values = TOC_conc
ds_all['TOC_c_heavy'].values = TOC_conc_heavy

# Get the surface layer
ds_all = ds_all.isel(siglay=0)

# # Get the surface layer
# ds_all = ds_all.isel(siglay=-1)
# print('mean depth', (ds_all['h']*ds_all['siglay']).mean().values)

# ds_all = ds_all.max(dim='siglay')

# Global min/max for consistent color scaling
# Not sure about compute() here
# vminmax = {
#     var: (ds_all[var].min().compute().item(), ds_all[var].max().compute().item())
#     for var in ['N_c', 'N_c_heavy', 'P_c', 'P_c_heavy', 'TOC_c', 'TOC_c_heavy']
# }
# # MAX
# vminmax = {}
# vminmax['N_c'] = (0, 0.4)
# vminmax['N_c_heavy'] = (0, 5)
# vminmax['P_c'] = (0, 0.1)
# vminmax['P_c_heavy'] = (0, 2)
# vminmax['TOC_c'] = (0, 1)
# vminmax['TOC_c_heavy'] = (0, 15)

# surface
vminmax = {}
vminmax['N_c'] = (0, 0.5)
vminmax['N_c_heavy'] = (0, 3)
vminmax['P_c'] = (0, 0.1)
vminmax['P_c_heavy'] = (0, 2)
vminmax['TOC_c'] = (0, 1)
vminmax['TOC_c_heavy'] = (0, 3)

# Setup figure and axes
fig, axs = plt.subplots(3, 2, figsize=(12, 12), sharex=True, sharey=True)
plots = {}
colorbars = {}

time_discharge = time[28] # 14th day
nt = len(time)

ds_all = ds_all.isel(time=slice(28, nt-1))  # Slice dataset to start from the 14th day
time = ds_all.time.values
nt = len(time)  # Get new number of timesteps after slicing


x, y = ds_all['x'], ds_all['y']
origin = [x[pipe_ind], y[pipe_ind]]

def init():
    ds = ds_all.isel(time=0)
    variables = ['N_c', 'N_c_heavy', 'P_c', 'P_c_heavy', 'TOC_c', 'TOC_c_heavy']
    titles = ['Nitrogen passive, g/m3', 'Nitrogen sinking, g/m3', 
              'Phosphorus passive, g/m3', 'Phosphorus sinking, g/m3', 'TOC passive, g/m3', 'TOC sinking, g/m3']
    cmaps = ['PuBu', 'PuBu', 'YlGn', 'YlGn', 'YlOrBr', 'YlOrBr']
    # cmaps = ['plasma', 'viridis', 'magma_r', 'magma_r']
    positions = [(0, 0), (0, 1), (1, 0), (1, 1), (2, 0), (2, 1)]

    for var, title, cmap, pos in zip(variables, titles, cmaps, positions):
        ax = axs[pos]

        plots[var] = ax.scatter(x, y, c=ds[var], cmap=cmap, s=1,
                                vmin=vminmax[var][0], vmax=vminmax[var][1])
        
        ax.set_title(title)
        ax.set_xlim(9.34e5, 9.365e5)
        ax.set_ylim(7.8535e6, 7.85525e6)
        ax.grid(True)
        ctx.add_basemap(ax, crs='EPSG:32633', source=ctx.providers.OpenStreetMap.Mapnik, attribution_size=6)

        # Add static colorbar
        colorbars[var] = fig.colorbar(plots[var], ax=ax, label=var)

    return plots.values()

def update(i):
    ds = ds_all.isel(time=i)
    current_time = time[i]
    current_date = fvcom_origin + timedelta(days=current_time)
    # time is in days since FVCOM origin
    delta_days = (current_time - time_discharge)

    for var in plots:
        plots[var].set_array(ds[var].values)

    if delta_days >= 0:
        fig.suptitle(f"Days since discharge: {delta_days:.2f} \n{current_date.strftime('%Y-%m-%d')}", fontsize=16)
    else:
        fig.suptitle(f"Days before discharge: {abs(delta_days):.2f} \n{current_date.strftime('%Y-%m-%d')}", fontsize=16)

    print(f"Frame {i} updated")
    return plots.values()

# Animate
anim = FuncAnimation(fig, update, init_func=init, frames=nt, blit=False)

# Save as MP4
writer = FFMpegWriter(fps=4, metadata=dict(artist='Python'), bitrate=1800)
anim.save("output_animation_surf3.mp4", writer=writer)

plt.close()
