import xarray as xr
import numpy as np
import matplotlib.pyplot as plt

ds = xr.open_dataset('forcing/adamselv_forcing_evap0.nc', decode_times=False)
# ds = xr.open_dataset('output/adamselv_forcing_TEST.nc', decode_times=False)


for k in ds.keys():
    print(k)

# # Convert time from Modified Julian Date (if needed)
# time_units = "days since 1858-11-17 00:00:00"  # Check metadata
# ds['time'] = xr.decode_cf(ds)['time']

print(ds)

def check_missing_values(ds):
    for var in ds.data_vars:
        missing_count = ds[var].isnull().sum().item()
        print(f"{var}: Missing values = {missing_count}")

# Run quality check
check_missing_values(ds)


variables_to_plot = ["U10", "V10", "air_temperature", "air_pressure", 
                     "relative_humidity", "short_wave", "long_wave", "precip"]

# fig, axes = plt.subplots(len(variables_to_plot), 1, figsize=(12, 15), sharex=True)

# for i, var in enumerate(variables_to_plot):
#     if var in ds:
#         if "node" in ds[var].dims:
#             ds[var].mean(dim="node").plot(ax=axes[i], label=var)
#         elif "nele" in ds[var].dims:
#             ds[var].mean(dim="nele").plot(ax=axes[i], label=var)
#         else:
#             axes[i].text(0.5, 0.5, f"{var} has no 'node' or 'nele' dimension", horizontalalignment='center')
#         axes[i].set_ylabel(var)
#         axes[i].legend()
#     else:
#         axes[i].text(0.5, 0.5, f"{var} not found", horizontalalignment='center')

# plt.xlabel("Time")
# plt.suptitle("Time Series of Atmospheric Forcing Variables")
# plt.tight_layout()
# plt.savefig("atmospheric_forcing_time_series.png", dpi=300, bbox_inches="tight")

if "SAT" in ds:
    plt.figure(figsize=(6, 10))
    sc = plt.scatter(ds.x, ds.y, c=ds.SAT.isel(time=0), cmap="coolwarm", s=1)
    plt.colorbar(sc, label="Surface Air Temperature (°C)")
    plt.title("Mean Surface Air Temperature")
    plt.savefig("atm_T_t0.png", dpi=300, bbox_inches="tight")

if "U10" in ds:
    plt.figure(figsize=(6, 10))
    sc = plt.scatter(ds.xc, ds.yc, c=ds.U10.isel(time=0), s=1)
    plt.colorbar(sc, label="U10 (m/s)")
    plt.title("U10")
    plt.savefig("atm_U10_t0.png", dpi=300, bbox_inches="tight")

if "precip" in ds:
    plt.figure(figsize=(6, 10))
    sc = plt.scatter(ds.x, ds.y, c=ds.precip.isel(time=0), s=1)
    plt.colorbar(sc, label="precipitation (m/s)")
    plt.title("precipitation")
    plt.savefig("atm_precip_t0.png", dpi=300, bbox_inches="tight")

if "evap" in ds:
    plt.figure(figsize=(6, 10))
    sc = plt.scatter(ds.x, ds.y, c=ds.evap.isel(time=0), s=1)
    plt.colorbar(sc, label="evaporation (m/s)")
    plt.title("evaporation")
    plt.savefig("atm_precip_evap.png", dpi=300, bbox_inches="tight")

