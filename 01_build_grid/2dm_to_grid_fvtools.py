import numpy as np
import sys
import os
import matplotlib.pyplot as plt

# THIS CREATES M.npy, WHICH IS USED TO GET M2.npy

sys.path.append('/Users/Admin/Documents/scripts/fvtools')

import fvtools.pre_pro.BuildCase as bc

# Set the current working directory to the script's directory
script_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(script_dir)

project = 'adamselv'

meshpath = '/Users/Admin/Documents/scripts/fvcom-work/Adamselv/adamselv_v01.2dm'
bathymetry = '/Users/Admin/Documents/scripts/fvcom-work/Adamselv/bathy_hybrid_adamselv_v01.npy'
sigma_file = '/Users/Admin/Documents/scripts/fvcom-work/Adamselv/input/adamselv_sigma.dat'
projection = 'epsg:32635' # UTM 35N
filename = '/Users/Admin/Documents/scripts/fvcom-work/Adamselv/input/adamselv_grd.dat'

ver_case = bc.main(meshpath, bathymetry,
                   casename='adamselv', 
                   sigma_file=sigma_file,
                   dm_projection = projection, 
                   depth_projection = projection, 
                   target_projection = 'epsg:32633',
                   rx0max = 0.18,
                   SmoothFactor = 0.18,)
# Potentially needs restarting with ver_case.main() to get past bad triangle

ver_case = ver_case.main()

# Write the bathymetry as for some reason 
M =  np.load('M.npy', allow_pickle=True).item()
dep = np.vstack([M['x'], M['y'], M['h']]).T
 
with open('input/adamselv_dep.dat', 'w') as f:
    for this_row in dep:
        f.write(f'{this_row[0]:06f} {this_row[1]:06f} {this_row[2]:06f}\n')

