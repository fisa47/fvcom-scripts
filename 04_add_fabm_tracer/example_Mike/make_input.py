import numpy as np
import datetime as dt
import netCDF4 as nc
from pyproj import Proj
import copy
import matplotlib.pyplot as plt

import PyFVCOM as pf
import river_funcs as rf


"""
Location
Adamselv: 29.5, 70.5 (SAMPLE)
"""

both_locations = {'hardanger1_v01': [29.5, 70.5]}
start_dates = {'exp1': dt.datetime(2013,8,1)}

river_basedir = '/home/michael/Projects/Vestlandpharma/preproc/flux'

release_length = dt.timedelta(minutes=15)
total_release = 7500
outflow_conc = 10
outflow_psu = 26

for hardanger_loc in [1,2]:
    for exp_name, start_date in start_dates.items():
        casename = f'hardanger{hardanger_loc}_v01'

        end_release = start_date + release_length
        end_run = dt.datetime(2020,1,1)

        outflow = total_release*(release_length/dt.timedelta(hours=1))  # Convert to m3/hr
        outflow_m3_per_sec = outflow/3600 # convert m3/hr to m3/sec

        grid_dir = '/home/michael/Projects/Vestlandpharma/preproc/grid'

        pipe_names = ['pipe1']
        #tri, nodes, x, y, z = pf.grid.read_fvcom_mesh(f'{grid_dir}/{casename}_grd.dat')

        out_file = nc.Dataset(f'/home/michael/nird/apn_backup/mib/Projects/Benchmark_hardanger/Benchmark_run_output/hardanger{hardanger_loc}/output/exp_1/hardanger{hardanger_loc}_v01_0001.nc')
        node_sig_deps = out_file['siglay_center'][:] * out_file['h_center'][:][np.newaxis,:]

        x = out_file['x'][:]
        y = out_file['y'][:]
        xc = out_file['xc'][:]
        yc = out_file['yc'][:]
        tri = np.asarray(out_file['nv'][:] - 1).T

        utm33 = Proj(proj = 'utm', zone = 33, ellps = 'WGS84', preserve_units = False)
        pipe_locs = np.asarray([utm33(this_loc[0], this_loc[1]) for this_loc in both_locations[casename]])

        pipe_eles = np.asarray([np.argmin(np.sqrt((x - pl[0])**2 + (y - pl[1])**2)) for pl in pipe_locs]) # Actually at the nodes
        pipe_deps = node_sig_deps[:,pipe_eles].T

        diameter = 3
        depth = 2
        outflow_deps = [depth - diameter/2, depth + diameter/2]

        pipe_dep_split = []

        for this_deps in pipe_deps:
            choose = np.logical_and(this_deps >= -outflow_deps[1], this_deps <= -outflow_deps[0])
            this_split = np.zeros(len(choose))
            this_split[choose] = 1/np.sum(choose)
            pipe_dep_split.append(this_split)

        pipe_dep_split = np.asarray(pipe_dep_split)
        pipe_dep_split_str = []

        for this_pipe in pipe_dep_split:
            this_str = ''
            for this_out in this_pipe:
                this_str = this_str + f'{this_out} '
            pipe_dep_split_str.append(this_str[0:-1])

        pipe_dt = [start_date + dt.timedelta(hours=int(i)) for i in np.arange(0, int((end_run - start_date).total_seconds()/(60*60))+1, 3)]

        # Make flux

        pipe_flux = np.asarray(np.ones([len(pipe_eles), len(pipe_dt)])*(outflow_m3_per_sec/len(pipe_eles))).T

        choose_stop = np.asarray(pipe_dt) >= end_release
        pipe_flux[choose_stop,:] = 0

        pipe_dt_hour = np.asarray([this_dt.hour for this_dt in pipe_dt])


        # Get salinity and temp from existing run by using the mean of the intake locations in the non-pipe run
        water_col_nc = nc.Dataset(f'water_column_hardanger{hardanger_loc}_pt1.nc')
        water_temp = np.mean(water_col_nc['temp'][:,0])
        water_salinity = np.mean(water_col_nc['salinity'][:,0])

        pipe_temp = np.ones(pipe_flux.shape)*water_temp
        pipe_salt = np.ones(pipe_flux.shape)*outflow_psu

        # Get the original river file to append these to

        original_file = f'{river_basedir}/hardanger{hardanger_loc}_river_fvtools.nc'
        original_nc = nc.Dataset(original_file)

        original_config = pf.read.get_river_config(f'{river_basedir}/hardanger{hardanger_loc}_v01_river_fvtools_node.nml')
        ref_date = dt.datetime(1858,11,17)
        original_dt = np.asarray([ref_date + dt.timedelta(days=float(t)) for t in original_nc['time'][:]])

        choose_original = np.logical_and(original_dt >= start_date, original_dt <= end_run)

        river_names_all = np.hstack([original_nc['river_names'][:], pipe_names]) 
        river_flux_all = np.hstack([original_nc['river_flux'][choose_original, :], pipe_flux])
        river_temp_all = np.hstack([original_nc['river_temp'][choose_original, :], pipe_temp])
        river_salt_all = np.hstack([original_nc['river_salt'][choose_original, :], pipe_salt])

        pipe_c1 = np.hstack([np.zeros(original_nc['river_flux'][choose_original, :].shape), np.ones((len(pipe_dt),1))*outflow_conc])

        if exp_name == 'baseline':
            pipe_c1[:] = 0
            pipe_flux[:] = 0

        # Add extra 2 timesteps for turning off the pipe

        add_dt = np.asarray([pipe_dt[0],pipe_dt[0] + release_length, pipe_dt[0] + release_length + dt.timedelta(seconds=15)]) 
        pipe_flux = np.hstack([outflow_m3_per_sec, outflow_m3_per_sec, np.zeros(len(pipe_dt))])

        pipe_dt = np.hstack([add_dt, pipe_dt[1:]])    
        river_flux_all = np.vstack([np.tile(river_flux_all[0,:], [3,1]), river_flux_all[1:,:]])
        river_temp_all = np.vstack([np.tile(river_temp_all[0,:], [3,1]), river_temp_all[1:,:]])
        river_salt_all = np.vstack([np.tile(river_salt_all[0,:], [3,1]), river_salt_all[1:,:]])
        river_flux_all[:,-1] = pipe_flux
        pipe_c1 = np.vstack([np.tile(pipe_c1[0,:], [3,1]), pipe_c1[1:,:]])

        river_locs_all = np.hstack([original_config['RIVER_GRID_LOCATION'][0:-1], pipe_eles])  
        river_vertical_dist_all = np.hstack([original_config['RIVER_VERTICAL_DISTRIBUTION'][0:-1], pipe_dep_split_str])

        # Write the river files
        tracer_dict = {}

        tracer_dict['river_names'] = river_names_all
        tracer_dict['no_rivers'] = river_flux_all.shape[1]  # 0 indexed
        tracer_dict['model_time'] = nc.date2num(pipe_dt, units='days since 1858-11-17 00:00:00')
        tracer_dict['flux'] = river_flux_all
        tracer_dict['temp'] = river_temp_all
        tracer_dict['salinity'] = river_salt_all

        tracer_dict['tracer1_c'] = pipe_c1

        rf.write_nc(tracer_dict, f'{casename}_{exp_name}_river_pipe.nc', 1)

        tracer_nml = {}
        tracer_nml['RIVER_GRID_LOCATION'] = river_locs_all
        tracer_nml['RIVER_NAME'] = tracer_dict['river_names']
        tracer_nml['RIVER_FILE'] = f'{casename}_river_pipe.nc'
        tracer_nml['RIVER_VERTICAL_DISTRIBUTION'] = river_vertical_dist_all

        rf.write_river_namelist(f'{casename}_{exp_name}_river_with_pipe.nml', tracer_nml)
         
        fig, ax = plt.subplots(figsize = [14,7])
        ax.plot(pipe_dt, pipe_c1[:,-1]*pipe_flux*1000)
        ax.set_ylabel('Pipe load (micro g)')
        ax.xaxis.set_major_locator(plt.MaxNLocator(6))
        fig.tight_layout()
        fig.savefig(f'pipe_outflow_load_{casename}_{exp_name}.png')
        plt.close()

