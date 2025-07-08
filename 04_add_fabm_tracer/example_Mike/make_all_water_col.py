import numpy as np
import netCDF4 as nc
import glob as gb
import xarray as xr
import datetime as dt
import time
import sys
import multiprocessing as mp
from pyproj import Proj

# Taken from PyFVCOM 

def proc_file(this_file_data):
    this_ind = this_file_data[0]
    this_file = this_file_data[1]
    this_nc = nc.Dataset(this_file)
    res_dict = {'id':this_ind}

    time_str_raw = this_nc.variables['Times'][:]
    res_dict['times'] = [dt.datetime.strptime(b''.join(this_str).decode('utf-8'), '%Y-%m-%dT%H:%M:%S.000000') for this_str in time_str_raw]

    node_vars_to_add = ['salinity', 'temp', 'zeta']
    for this_var in node_vars_to_add:
        res_dict[this_var] = this_nc[this_var][...,chosen_node]
        
    ele_vars_to_add = ['u', 'v']
    for this_var in ele_vars_to_add:
        res_dict[this_var] = this_nc[this_var][...,chosen_ele]

    return res_dict

def array_from_dict(in_d):
    d_keys = list(in_d.keys())
    d_keys.sort()
    out = []
    for this_key in d_keys:
        out.append(in_d[this_key])
    out = np.vstack(out)
    return np.squeeze(out)

if __name__ == "__main__":
    basedir = '/nird/projects/NS9067K/apn_backup/FVCOM/Peygham/Svalbard/Svalbard_newgrid'

    met_file = '/nird/projects/NS9067K/apn_backup/FVCOM/Peygham/Svalbard/Svalbard3D/Sval3/input/Sval3_wnd.nc'
    nprocs = 8

    all_locs = {'awake': [16.1747, 76.984], 'PEV':[11.919042, 78.930447],
                'FS':[8.8627, 78.8328], 'PI1':[13.3866, 78.1816], 'PI2':[13.5231, 78.0609],
                'PI3':[14.4172, 78.1383], 'PI4':[15.3312, 78.2474], 'PI5':[17.3538, 78.4433],
                'Isfjorden':[13.39,78.18],  'S1':[13.9484, 76.4381], 'KF':[11.8238, 78.9589], 'VF':[16.61256, 77.82918]}

    for loc_name, loc in all_locs.items():
        print(loc_name, flush=True)
        lon = loc[0]
        lat = loc[1]
        test_file = '/nird/projects/NS9067K/apn_backup/FVCOM/Peygham/Svalbard/Svalbard_newgrid/2016_noOBSICE/output201601/Sval3_0001.nc'
        test_nc = nc.Dataset(test_file)

        utm33 = Proj(proj = 'utm', zone = 33, ellps = 'WGS84', preserve_units = False)
        x, y = utm33(lon, lat)

        chosen_node = np.argmin((test_nc['x'][:] - x)**2 + (test_nc['y'][:] - y)**2)
        chosen_ele = np.argmin((test_nc['xc'][:] - x)**2 + (test_nc['yc'][:] - y)**2)

        vars_to_add = ['times', 'salinity', 'temp', 'zeta', 'u', 'v']
        out_dict = {}
        for this_var in vars_to_add:
            out_dict[this_var] = []

        for this_year in [2016,2017]:
            print(this_year, flush=True)
            for this_month in np.arange(1,13):
                print(this_month, flush=True)
                this_monthdir = f'{basedir}/{this_year}_noOBSICE/output{this_year}{this_month:02d}'
            
                raw_flist = gb.glob(f'{this_monthdir}/Sval3_00*.nc')
                raw_flist.sort()
                filelist = [(i, filename) for i, filename in enumerate(raw_flist)]

                pool = mp.Pool(nprocs)
                res = pool.map(proc_file, filelist)

                sort_dict = {}
                # This next code is ineffecient, I'm sure theres a better way to do it...
                for this_var in vars_to_add:
                    sort_dict[this_var] = {}
                
                for this_row in res:
                    for this_var in vars_to_add:
                        sort_dict[this_var][this_row['id']] = this_row[this_var]
                 
                for this_var in vars_to_add:
                    out_dict[this_var].append(array_from_dict(sort_dict[this_var]))

        #np.save('out_dict.npy', out_dict)

        # Trim off overlapping days
        new_out_dict = {}
        for var_name, this_var in out_dict.items():
            new_var = []
            if this_var[0].shape[-1] == 34:    
                for this_entry in this_var:
                    new_var.append(this_entry[0:-24,:])
                new_out_dict[var_name] = np.vstack(new_var)
            elif this_var[0].shape[-1] == 24:
                for this_entry in this_var:
                    new_var.append(this_entry[0:-1,:].ravel())
                new_out_dict[var_name] = np.hstack(new_var)
            else:
                raise ValueError

        ds = xr.Dataset(data_vars={'salinity':(['time', 'depth'], new_out_dict['salinity']), 
            'temp':(['time', 'depth'], new_out_dict['temp']),
               'u':(['time', 'depth'], new_out_dict['u']), 
               'v':(['time', 'depth'], new_out_dict['v']), 
               'zeta':(['time'], new_out_dict['zeta']),
               'siglay':(['depth'], np.squeeze(test_nc['siglay'][:,chosen_node])),
                    }, coords={'depth':np.arange(0, new_out_dict['salinity'].shape[1]), 'time':new_out_dict['times'], 'h':test_nc['h'][chosen_node]})

        ds.time.encoding['units'] = 'seconds since 2016-01-01'
        ds.to_netcdf(f'{loc_name}_water_column.nc')

