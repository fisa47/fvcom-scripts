import xarray as xr
import matplotlib.pyplot as plt
import contextily as ctx
import numpy as np
from datetime import datetime, timedelta

fvcom_origin = datetime(1858, 11, 17)

pipe_discharge = 5.3 / 60
yr_vol = pipe_discharge * 86400 * 365
to_ug_per_L = 1e6 / yr_vol
TOC_ug_L_part = 338557 * to_ug_per_L
pipe_ind = 54468

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

# Load dataset
ds = xr.open_mfdataset('../run-output/pipe_ave/adamselv_v01_000*.nc', decode_times=False, concat_dim='time', combine='nested')

xpipe = ds['x'].isel(node=pipe_ind).compute().item()
ypipe = ds['y'].isel(node=pipe_ind).compute().item()

# Extract variables
area = ds['art1'].values
depth = ds['h'].values
surface_elevation = ds['zeta'].values
thickness = ds['siglev'].values
time = ds.time.values
tri = ds['nv'].isel(time=0).values.T - 1
vol = unstructured_grid_volume(area, depth, surface_elevation, thickness)

# Compute corrected TOC
TOC_conc_part = correct_tracer_concentration(ds, 'tracer2_c', vol, time)
ds['TOC_c_part'] = ds['tracer2_c'] * 0
ds['TOC_c_part'].values = TOC_conc_part

target_ind = np.array([53601, 53763, 53918, 54066, 54208, 54402, 54468, 54531, 54591, 54694, 54744])
dsout = ds.sel(node=target_ind)

toc = dsout['TOC_c_part']
x = dsout.x.values - xpipe
y = dsout.y.values - ypipe
yout = np.tile(y[:, np.newaxis], 30).T
zout = dsout['siglay'].values * dsout['h'].isel(time=0).values

# Plotting setup
times_to_plot = [50, 100, 150, 200, 250]
titles = []
for t in times_to_plot:
    current_time = time[t]
    current_date = fvcom_origin + timedelta(days=float(current_time))
    days_since_discharge = float(current_time - time[28])
    titles.append(f"TOC, μg/L (Day {days_since_discharge:.0f})")

fig, axs = plt.subplots(2, 3, figsize=(14, 6))

# Map
ax = axs[0, 0]
ax.set_xlim(9.352e5, 9.3625e5)
ax.set_ylim(7.854e6, 7.8546e6)
ax.scatter(x + xpipe, y + ypipe, label='transect nodes')
ax.scatter(xpipe, ypipe, label='source')
ax.legend()
ctx.add_basemap(ax, crs='EPSG:32633', source=ctx.providers.CartoDB.Positron, attribution_size=6)
ax.set_title("Transect map")

# Plot TOC profiles
for i, t in enumerate(times_to_plot):
    row = (i + 1) // 3
    col = (i + 1) % 3
    ax_t = axs[row, col]
    toc_plot = ax_t.tricontourf(yout.ravel(), zout.ravel(), toc.isel(time=t).values.ravel(), levels=50, cmap='YlOrBr')
    ax_t.set_title(titles[i])
    ax_t.set_xlabel('Distance from pipe (m)')
    ax_t.set_ylabel('Depth (m)')
    fig.colorbar(toc_plot, ax=ax_t, fraction=0.046, pad=0.04)

plt.tight_layout()
plt.savefig('transect_multi_time.png', dpi=200, bbox_inches='tight')
