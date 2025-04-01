import xarray as xr
import numpy as np
import matplotlib.pyplot as plt

"""
This script modifies an existing NetCDF file by setting all evaporation values to zero and adjusting the time values.

Steps:
1. Open the NetCDF dataset using xarray.
2. Set all values of the 'evap' variable to zero.
3. Adjust the 'time' variable by subtracting 1826 days (5 year from 2013 to 2018).
4. Save the modified dataset to a new NetCDF file.

Example:
python fix_evaporation.py

Dependencies:
- xarray
- numpy
- matplotlib
"""

ds = xr.open_dataset('output/adamselv_forcing_201802-06.nc', decode_times=False)
shape = ds['evap'].shape
ds['evap'].values = np.zeros(shape)  # set all values to zero

ds = ds.assign_coords(time=ds['time'].values - 1826)  # subtract 5 years
ds['Itime'].values = ds['Itime'].values - 1826  # substract 5 years
ds.to_netcdf('adamselv_forcing_feb-jun_evap0.nc')