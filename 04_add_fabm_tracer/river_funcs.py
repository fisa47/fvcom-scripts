import numpy as np
import netCDF4
import os
from datetime import datetime
import PyFVCOM as pf

def write_nc(data, out_name, no_tracers):
    """
    Write the riverdata.nc file for river forcing
    """
    # Initialize file
    d = netCDF4.Dataset(out_name, 'w')

    # Set dimensions
    d.createDimension('time', None)
    d.createDimension('rivers', data['no_rivers'])
    d.createDimension('DateStrLen', 26)
    d.createDimension('namelen', 80)

    # Add netcdf information
    d.source      = 'fvtools functions used seperately'
    d.history     = 'Created '+ datetime.now().strftime('%Y-%m-%d at %H:%M h')+' by '+os.getlogin()
    d.description = 'River forcing (temperature and runoff) for FVCOM 4.x'

    # Create variables:
    # - time
    time = d.createVariable('time', 'single', ('time',))
    time.long_name   = 'time'
    time.units       = 'days since '+str(datetime(1858, 11, 17, 0, 0, 0))
    time.format      = 'modified julian day (MJD)'
    time.time_zone   = 'UTC'

    # - Itime
    Itime = d.createVariable('Itime', 'int32', ('time',))
    Itime.long_name   = 'integer days'
    Itime.units       = 'days since '+str(datetime(1858, 11, 17, 0, 0, 0))
    Itime.format      = 'modified julian day (MJD)'
    Itime.time_zone   = 'UTC'

    # - Itime2
    Itime2           = d.createVariable('Itime2', 'int32', ('time',))
    Itime2.long_name = 'integer milliseconds'
    Itime2.units     = 'msec since 00:00:00'
    Itime2.time_zone = 'UTC'

    # - river_flux
    flux             = d.createVariable('river_flux', 'single', ('time','rivers'))
    flux.long_name   = 'river runoff volume flux, m**-3 s**-1'
    flux.units       = 'm^3s^-1'

    # - river_temp
    temp             = d.createVariable('river_temp', 'single', ('time','rivers'))
    temp.long_name   = 'river runoff temperature'
    temp.units       = 'Celsius'

    # - river_salt
    salt             = d.createVariable('river_salt', 'single', ('time','rivers'))
    salt.long_name   = 'river runoff salinity'
    salt.units       = 'PSU'

    # - river_names
    names            = d.createVariable('river_names', 'S1', ('rivers', 'namelen'))

    # - river_flux
    flux_list = []
    for i in np.arange(1,no_tracers+1):
        this_flux = d.createVariable(f'tracer{i}_c', 'single', ('time','rivers'))
        this_flux.long_name   = f'tracer{i} concentration'
        this_flux.units       = 'kg/m^3'

        flux_list.append(this_flux)

    # Dump data:
    salt[:]   = data['salinity']
    temp[:]   = data['temp']
    flux[:]   = data['flux']
    time[:]   = data['model_time']
    Itime[:]  = data['model_time']
    Itime2[:] = np.round((data['model_time'] - np.floor(data['model_time'])) * 60 * 60 * 1000, decimals = 0)*24

    for i, this_flux in enumerate(flux_list):
        this_flux[:] = data[f'tracer{i+1}_c']

    # Dump river names
    # --> Make sure that each rivername has 80 character
    _names = []
    names._Encoding = 'ascii'
    for i, name in enumerate(data['river_names']):
        if len(name) > 80:
            this_name = name[:80]
        else:
            this_name  = name + (80-len(name))*' '
        names[i,:] = np.array(this_name, dtype = 'S80')

    d.close()

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
                                  pf.preproc.NameListEntry('RIVER_FILE', conf_dict['RIVER_FILE']),
                                  pf.preproc.NameListEntry('RIVER_GRID_LOCATION', conf_dict['RIVER_GRID_LOCATION'][ri] + 1, 'd'),
                                  pf.preproc.NameListEntry('RIVER_VERTICAL_DISTRIBUTION', conf_dict['RIVER_VERTICAL_DISTRIBUTION'][ri])]}
        pf.preproc.write_model_namelist(output_file, namelist, mode='a')

