import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegWriter
import contextily as ctx
import pandas as pd
import numpy as np
from datetime import datetime, timedelta
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D

# FVCOM origin date
fvcom_origin = datetime(1858, 11, 17)

pipe_discharge = 5.3/60  # m3/s, from Grieg email

# Nutrients data from
# Notat 62802, Sak: Søknad om endring i utslippstillatelse Adamselv Settefisk
# Table 5
yr_vol = pipe_discharge * 86400 * 365  # m3
to_ug_per_L   = 1e6 / yr_vol # conversion factor

# NEW NUMBERS # µg/L
# Grieg Adamselv - Utslippsberegning fra Skretting_TIRO #4
# https://akvaplan.sharepoint.com/:x:/r/sites/projects3/6407/_layouts/15/Doc.aspx?sourcedoc=%7B851DCCA6-5FD9-407B-BB10-97A5BD3436F6%7D&file=Grieg%20Adamselv%20-%20Utslippsberegning%20fra%20Skretting_TIRO%20%234.xlsx&nav=MTBfezBCMTk0OTRDLTlCMUQtNERBMC1CMzZDLTZEMjRENUVGODQ3Nn1fezBDQ0REOTk3LTQ5M0QtNEQxRi04OEY4LTk2ODI0RkMwRTM1MX0&action=default&mobileredirect=true
N_ug_L_diss = 98595 * to_ug_per_L  
P_ug_L_diss = 4347 * to_ug_per_L  

N_ug_L_part = 36762 * to_ug_per_L  
P_ug_L_part = 22533 * to_ug_per_L  
TOC_ug_L_part = 338557 * to_ug_per_L  


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
    print(f"Tracer {tracer_name} corrected")
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
print(tri.shape)

# Calculate volume and mass
vol = unstructured_grid_volume(area, depth, surface_elevation, thickness)

N_conc_diss = correct_tracer_concentration(ds_all, 'tracer1_c', vol, time, 150, substance_concentration=N_ug_L_diss)
N_conc_part = correct_tracer_concentration(ds_all, 'tracer2_c', vol, time, 150, substance_concentration=N_ug_L_part)

P_conc_diss = correct_tracer_concentration(ds_all, 'tracer1_c', vol, time, 10, substance_concentration=P_ug_L_diss)
P_conc_part = correct_tracer_concentration(ds_all, 'tracer2_c', vol, time, 10, substance_concentration=P_ug_L_part)

# there is no dissolved TOC
TOC_conc_part = correct_tracer_concentration(ds_all, 'tracer2_c', vol, time, 1500, substance_concentration=TOC_ug_L_part)

# Create new variables for corrected concentrations
ds_all['N_c_diss'] = ds_all['tracer1_c'] * 0
ds_all['N_c_part'] = ds_all['tracer1_c'] * 0

ds_all['P_c_diss'] = ds_all['tracer1_c'] * 0
ds_all['P_c_part'] = ds_all['tracer1_c'] * 0

ds_all['TOC_c_part'] = ds_all['tracer1_c'] * 0

# Assign corrected concentrations to new variables
ds_all['N_c_diss'].values = N_conc_diss
ds_all['N_c_part'].values = N_conc_part

ds_all['P_c_diss'].values = P_conc_diss
ds_all['P_c_part'].values = P_conc_part

ds_all['TOC_c_part'].values = TOC_conc_part

# Get the max
# ds_all = ds_all.max(dim='siglay', skipna=True).max(dim='time')
ds_all = ds_all.isel(siglay=0).mean(dim='time')


# surface
vminmax = {}
vminmax['N_c_diss'] = (150, N_ug_L_diss/500)
vminmax['N_c_part'] = (150, N_ug_L_part/20)
vminmax['P_c_diss'] = (0, 8)
vminmax['P_c_part'] = (10, P_ug_L_part/20)

vminmax['TOC_c_part'] = (1500, TOC_ug_L_part/20)

# https://www.statsforvalteren.no/siteassets/fm-agder/dokument-agder/miljo-og-klima/mudring-utfylling-og-dumping/2018---02-veileder-klassifisering-av-miljotilstand-i-vann-27.10.20-1.pdf
# p. 173, Table 9.26, S=18, Sommer (Juni, Juli, Aug)
levels = {}
                  #  Svært god, God, Moderat, Dårlig
levels['N_c_diss'] = [12, 23, 65, 250]  # Nitrat + nitritt
levels['N_c_part'] = [12, 23, 65, 250]
levels['P_c_diss'] = [3.5, 7, 16, 50]
levels['P_c_part'] = [3.5, 7, 16, 50]

levels['TOC_c_part'] = [1500, 4300, 1e4, 1e5]  # just typical values

# Setup figure and axes
fig, axs = plt.subplots(3, 2, figsize=(12, 12), sharex=True, sharey=True)
plots = {}
colorbars = {}

# time_discharge = time[28] # 14th day
# nt = len(time)

# ds_all = ds_all.isel(time=slice(28, nt-1))  # Slice dataset to start from the 14th day
# time = ds_all.time.values
# nt = len(time)  # Get new number of timesteps after slicing

# # Get timestep index where N_c is maximal over all elements and time
# max_light = ds_all['N_c_diss'].max(dim='node').argmax(dim='time').compute().item()
# max_part = ds_all['N_c_part'].max(dim='node').argmax(dim='time').compute().item()

# snapshot_steps = [0, max_light, max_part, nt-1]  # First, max, last timestep
# fnames = ['start', 'max_light', 'max_part', 'end']


x, y = ds_all['x'], ds_all['y']
origin = [x[pipe_ind], y[pipe_ind]]
extent = 1500

fnames = ['BOTTOM']

# for i, fname in zip(snapshot_steps, fnames):
for fname in fnames:
    print(f"Plotting {fname}")
    fig, axs = plt.subplots(3, 2, figsize=(12, 12), sharex=True, sharey=True)
    ds = ds_all # .isel(time=i)
    # current_time = time[i]
    # current_date = fvcom_origin + timedelta(days=float(current_time))
    # delta_days = current_time - time_discharge

    variables = ['N_c_diss', 'N_c_part', 'P_c_diss', 'P_c_part', 'TOC_c_part']
    titles = ['Dissolved Nitrogen passive, μg/L', 'Particulate Nitrogen sinking, μg/L',
              'Dissolved Phosphorus, μg/L', 'Paticulate Phosphorus, μg/L', ' Particulate TOC, μg/L']
    
    # This does not work
    cmaps = [plt.colormaps["PuBu"].with_extremes(under=None),
             plt.colormaps["PuBu"].with_extremes(under=None),
             plt.colormaps["YlGn"].with_extremes(under=None),
             plt.colormaps["YlGn"].with_extremes(under=None),
             plt.colormaps["YlOrBr"].with_extremes(under=None),
             ]
    positions = [(0, 0), (0, 1), (1, 0), (1, 1), (2, 1)]

    for var, title, cmap, pos in zip(variables, titles, cmaps, positions):
        ax = axs[pos]
        # sc = ax.scatter(x, y, c=ds[var].values, cmap=cmap, s=1,)
        
        sc = ax.tripcolor(x, y, tri, ds[var].values,
                          shading='flat',
                          cmap=cmap,
                          edgecolors='none',
                          vmin=vminmax[var][0],
                          vmax=vminmax[var][1])

        if var != 'TOC_c_part':
            contour = ax.tricontour(x, y, tri, ds[var].values,
                                levels=levels[var],
                                linewidths=2, colors=['dodgerblue', 'green', 'gold', 'orange'], alpha=1)
        ax.set_title(title)
        ax.set_xlim(9.34e5, 9.365e5)
        ax.set_ylim(7.8535e6, 7.85525e6)
        ax.grid(True)
        ctx.add_basemap(ax, crs='EPSG:32633', source=ctx.providers.OpenStreetMap.Mapnik, attribution_size=6)

        if var not in ['N_c_part', 'P_c_part', 'TOC_c_part']:
            # legend only for dissolved N and P
            legend_elements = [
                Line2D([0], [0], color='dodgerblue', lw=2, label= f'Svært god (< {levels[var][0]} μg/L)'),
                Line2D([0], [0], color='green', lw=2,      label=f'God ({levels[var][0]} - {levels[var][1]} μg/L)'),
                Line2D([0], [0], color='gold', lw=2,       label=f'Moderat ({levels[var][1]} - {levels[var][2]} μg/L)'),
                Line2D([0], [0], color='orange', lw=2,     label=f'Dårlig {levels[var][2]} - {levels[var][3]} μg/L)'),
            ]

            # Add legend to the lower-left corner of the figure
            ax.legend(handles=legend_elements, loc='upper left', frameon=True, fontsize=12,
                    title_fontsize=14, facecolor='0.95', edgecolor='0.5')

        fig.colorbar(sc, ax=ax, label=var)

    # goodbye axis
    axs[2, 0].remove()

    # # Set a time-aware supertitle
    # if delta_days >= 0:
    #     fig.suptitle(f"Days since discharge: {delta_days:.2f} \n{current_date.strftime('%Y-%m-%d')}", fontsize=16)
    # else:
    #     fig.suptitle(f"Days before discharge: {abs(delta_days):.2f} \n{current_date.strftime('%Y-%m-%d')}", fontsize=16)


    # Save to file
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(f"snapshot_AVE_surf.png", dpi=200)
    plt.close()
