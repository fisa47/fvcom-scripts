import xarray as xr
import matplotlib.pyplot as plt
import contextily as ctx
import numpy as np
from datetime import datetime, timedelta
from matplotlib.lines import Line2D

# Constants
fvcom_origin = datetime(1858, 11, 17)
pipe_discharge = 5.3 / 60
yr_vol = pipe_discharge * 86400 * 365
to_ug_per_L = 1e6 / yr_vol
TOC_ug_L_part = 338557 * to_ug_per_L
pipe_ind = 54469

# Load dataset
ds_all = xr.open_mfdataset('../run-output/pipe_ave/adamselv_v01_000*.nc', decode_times=False, concat_dim='time', combine='nested')

# Volume calculation
def unstructured_grid_volume(area, depth, surface_elevation, thickness):
    dz = np.abs(np.diff(thickness, axis=0))
    volume = (area * (surface_elevation + depth))
    depth_volume = volume[:, np.newaxis, :] * dz[np.newaxis, ...]
    return depth_volume

# Correction function
def correct_tracer_concentration(ds, tracer_name, vol, time,
                                 pipe_concentration=1500, pipe_discharge=5.3/60,
                                 substance_concentration=TOC_ug_L_part):
    tracer_mass = vol * ds[tracer_name].values
    tracer_mass_discharged = pipe_concentration * pipe_discharge * (np.abs(time[0] - time) * 86400)
    int_mass = tracer_mass.sum(axis=(-1, -2))
    calibration = tracer_mass_discharged / int_mass
    calibration[0] = 1
    corrected_mass = tracer_mass * calibration[:, np.newaxis, np.newaxis]
    substance_mass_discharged = substance_concentration * pipe_discharge * (np.abs(time[0] - time) * 86400)
    int_corrected_mass = corrected_mass.sum(axis=(-1, -2))
    calibration_substance = substance_mass_discharged / int_corrected_mass
    corrected_substance_mass = corrected_mass * calibration_substance[:, np.newaxis, np.newaxis]
    tracer_conc = corrected_substance_mass / vol
    print(f"Tracer {tracer_name} corrected")
    return tracer_conc

# Extract variables
area = ds_all['art1'].values
depth = ds_all['h'].values
surface_elevation = ds_all['zeta'].values
thickness = ds_all['siglev'].values
time = ds_all.time.values
tri = ds_all['nv'].isel(time=0).values.T - 1
vol = unstructured_grid_volume(area, depth, surface_elevation, thickness)

# Compute corrected TOC
TOC_conc_part = correct_tracer_concentration(ds_all, 'tracer2_c', vol, time)
ds_all['TOC_c_part'] = ds_all['tracer2_c'] * 0
ds_all['TOC_c_part'].values = TOC_conc_part

# Setup plot parameters
x, y = ds_all['x'], ds_all['y']
vmin, vmax = 1500, TOC_ug_L_part / 5
levels = [1500, 4300, 10000, 100000]

# Choose timesteps
print("Available timesteps:", len(time))
timesteps = [59, 89, 119, 149, 
             179, 209, 249, len(time) - 1]
titles = []
for t in timesteps:
    current_time = time[t]
    current_date = fvcom_origin + timedelta(days=float(current_time))
    days_since_discharge = float(current_time - time[28])
    titles.append(f"TOC, μg/L (Day {days_since_discharge:.0f})")

# Create figure
fig, axs = plt.subplots(4, 2, figsize=(10, 15), sharex=True, sharey=True)

for ax, t, title in zip(axs.flat, timesteps, titles):
    ds_t = ds_all.isel(siglay=0, time=t)
    sc = ax.tripcolor(x, y, tri, ds_t['TOC_c_part'].values,
                      shading='flat', cmap='YlOrBr', vmin=vmin, vmax=vmax)
    ax.set_title(title)
    ax.set_xlim(9.35e5, 9.365e5)
    ax.set_ylim(7.8540e6, 7.8550e6)
    ax.grid(True)
    ctx.add_basemap(ax, crs='EPSG:32633', source=ctx.providers.OpenStreetMap.Mapnik)
    fig.colorbar(sc, ax=ax, label='TOC, μg/L')

fig.tight_layout(rect=[0, 0.05, 1, 1])
fig.savefig("TOC_surface_snapshots.png", dpi=200)
plt.close()
