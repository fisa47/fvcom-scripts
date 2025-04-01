import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegWriter

# Load dataset
ds_all = xr.open_dataset('../run-output/pipe_ave/adamselv_v01_0001.nc', decode_times=False).isel(siglay=1)
print(ds_all)

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
