import numpy as np
import sys
import os
import matplotlib.pyplot as plt

sys.path.append('/Users/Admin/Documents/scripts/fvtools')

import fvtools.pre_pro.BuildCase as bc

# Set the current working directory to the script's directory
script_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(script_dir)

import fvtools.pre_pro.BuildRivers as br
br.main('2013-08-01-00', '2013-10-30-00', [227, 228, 229, 230], mesh_dict = 'M2.npy', temp='Finnmark_temperatures.npy')

