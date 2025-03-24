import xarray as xr
import numpy as np
import matplotlib.pyplot as plt

ds = xr.open_dataset('forcing/adamselv_forcing_2013.nc', decode_times=False)
# ds = xr.open_dataset('output/adamselv_forcing_TEST.nc', decode_times=False)
shape = ds['evap'].shape
ds['evap'].values = np.zeros(shape) # set all values to zero
ds.to_netcdf('forcing/adamselv_forcing_evap0.nc')