import numpy as np
import netCDF4 as nc
import sys

sys.path.append('/Users/Admin/Documents/scripts/pyfvcom')
import PyFVCOM as pf

casename = f'adamselv_v01'
out_file = nc.Dataset(f'adamselv_v01_0001.nc')
node_sig_deps = out_file['siglay_center'][:] * out_file['h_center'][:][np.newaxis,:]

x = out_file['x'][:]
y = out_file['y'][:]
xc = out_file['xc'][:]
yc = out_file['yc'][:]
tri = np.asarray(out_file['nv'][:] - 1).T

original_config = pf.read.get_river_config(f'input/RiverNamelist.nml')
old_locs = original_config['RIVER_GRID_LOCATION']

ob = np.hstack(pf.grid.get_boundary_polygons(tri))

new_locs = []

for this_loc in old_locs[0:-1]:
	new_locs.append(ob[np.argmin(np.sqrt((x[ob] - xc[this_loc])**2 + (y[ob] - yc[this_loc])**2))])

new_locs.append(np.argmin(np.sqrt((x - xc[old_locs[-1]])**2 + (y - yc[old_locs[-1]])**2)))
new_locs = np.asarray(new_locs)


# Get rid of overlaps



while len(new_locs) != len(np.unique(new_locs)):
	u,c = np.unique(new_locs, return_counts=True)

	to_sort = u[c>1]

	for this_node in to_sort:
		replace = np.squeeze(np.where(new_locs == this_node))[1]
		this_loc = new_locs[replace]
		dists = np.sqrt((x[ob] - xc[this_loc])**2 + (y[ob] - yc[this_loc])**2)
		dists[dists == np.min(dists)] = 10*10**10
		new_locs[replace] = ob[np.argmin(dists)]



def write_river_namelist(output_file, conf_dict):
    """
    Write an FVCOM river namelist file.

    Parameters
    ----------
    output_file : str, pathlib.Path
        Output file to which to write the river configuration.
    forcing_file : str, pathlib.Path
        File from which FVCOM will read the river forcing data.
    vertical_distribution : str, optional
        Vertical distribution of river input. Defaults to 'uniform'.

    """
    for ri in np.arange(0,len(conf_dict['RIVER_GRID_LOCATION'])):
        namelist = {'NML_RIVER': [pf.preproc.NameListEntry('RIVER_NAME', conf_dict['RIVER_NAME'][ri]),
                                  pf.preproc.NameListEntry('RIVER_FILE', conf_dict['RIVER_FILE'][ri]),
                                  pf.preproc.NameListEntry('RIVER_GRID_LOCATION', conf_dict['RIVER_GRID_LOCATION'][ri] + 1, 'd'),
                                  pf.preproc.NameListEntry('RIVER_VERTICAL_DISTRIBUTION', conf_dict['RIVER_VERTICAL_DISTRIBUTION'][ri])]}
        pf.preproc.write_model_namelist(output_file, namelist, mode='a')

original_config['RIVER_GRID_LOCATION'] = new_locs

write_river_namelist(f'input/{casename}_river_node.nml', original_config)

