import xarray as xr
import numpy as np
import matplotlib.pyplot as plt

data_url = "https://thredds.met.no/thredds/dodsC/meps25epsarchive/2018/08/01/meps_mbr0_extracted_backup_2_5km_20180801T00Z.nc"

# Open dataset
try:
    ds = xr.open_dataset(data_url)
    print("Dataset loaded successfully.")
except Exception as e:
    print(f"Failed to load dataset: {e}")
    exit()

print("Available variables:", ds.data_vars.keys())

def check_missing_values(ds):
    for var in ds.data_vars:
        missing_count = ds[var].isnull().sum().item()
        print(f"{var}: Missing values = {missing_count}")

check_missing_values(ds)

# Variable to plot
var_name = "water_evaporation_amount"

if var_name in ds:
    evaporation = ds[var_name]
    
    # Time-averaged
    plt.figure(figsize=(8, 6))
    evaporation.mean(dim="time").plot()
    plt.title("Mean Water Evaporation Amount")
    plt.colorbar(label="Evaporation (kg/m²)")
    plt.savefig("evaporation.png", dpi=300, bbox_inches="tight")
    plt.show()
else:
    print(f"Variable {var_name} not found in dataset.")