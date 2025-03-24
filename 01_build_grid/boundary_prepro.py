
import numpy as np
import sys
import os
import matplotlib

matplotlib.use('TkAgg')
sys.path.append('/Users/Admin/Documents/scripts/pyfvcom')

import matplotlib.pyplot as plt

import PyFVCOM as pf
from pyproj import Proj, transform

# Set the current working directory to the script's directory
script_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(script_dir)

# Tidy up open boundary

parent_grid = 'finntest/finntest_grd.dat'

utm33 = Proj(proj = 'utm', zone = 33, ellps = 'WGS84', preserve_units = True)
wgs84 = Proj(proj='latlong', datum='WGS84')

parent_file = f'./{parent_grid}'
tri, nodes, x, y, z = pf.grid.read_fvcom_mesh(parent_file)

x, y = transform(wgs84, utm33, x, y)

M =  np.load('M.npy', allow_pickle=True).item()

plt.triplot(x, y, tri, c='lightgray', linewidth=0.2)
plt.triplot(M['x'], M['y'], M['nv'], c='darkgray', linewidth=0.2)
plt.scatter(M['x'][M['obc_nodes']], M['y'][M['obc_nodes']], c='darkgreen', s=1, zorder=2)
plt.savefig('obc_nodes.png', dpi=300)

outer_bound = pf.grid.get_boundary_polygons(M['nv'])[0]

[70.9506,26.6673]
[70.9447,27.2379]

ns_dict_xy = {1: [
                  [943253.4,7916180.2],
                  [922692.8, 7912738.7],
                  ]}

"""
# ooops wrote the wrong way round
ns_dict_xy = {}
for this_key, this_val in ns_dict_xy_wrong.items():
    ns_dict_xy[this_key] = [this_val[1], this_val[0]]
"""

ns_dict = {}
for this_ns, this_xy in ns_dict_xy.items():
    end_pts = []
    for this_pt in this_xy:
        end_pts.append(np.argmin(np.sqrt((M['x'] - this_pt[0])**2 + (M['y'] - this_pt[1])**2)))

    end_ob = []
    for this_ind in end_pts:
        end_ob.append(np.where(outer_bound == this_ind)[0][0])

    end_ob[1] = end_ob[1] + 1
    if end_ob[1] < end_ob[0]:
        ns_dict[this_ns] = np.hstack([outer_bound[end_ob[0]:], outer_bound[0:end_ob[1]]])
    else:
        ns_dict[this_ns] = outer_bound[end_ob[0]:end_ob[1]]

plt.figure(figsize=(10, 3))
plt.triplot(x,y,tri,c='lightgray')
plt.triplot(M['x'],M['y'],M['nv'],c='darkgray')
for this_b in ns_dict.values():
    plt.scatter(M['x'][this_b], M['y'][this_b], c='tomato', s=2, zorder=2)

# plt.xlim(487000, 509000)
# plt.ylim(7870000, 7873000)
plt.savefig('obc_nodes_fixed.png', dpi=300)

# And re-write files
M['nodestrings'] = [this_str for this_str in ns_dict.values()]
M['obc_nodes'] = np.hstack(M['nodestrings'])
np.save('M2.npy', M)

def write_obc_file(obc_nodes, obc_types, obc_file):
    number_of_nodes = len(obc_nodes)

    with open(str(obc_file), 'w') as f:
        f.write('OBC Node Number = {:d}\n'.format(number_of_nodes))
        for count, node, obc_type in zip(np.arange(number_of_nodes) + 1, obc_nodes, obc_types):
            f.write('{} {:d} {:d}\n'.format(count, int(node) + 1, int(obc_type)))

write_obc_file(M['obc_nodes'], np.ones(len(M['obc_nodes'])), './input/fvtools_obc.dat')