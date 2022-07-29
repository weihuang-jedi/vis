#=========================================================================
import sys
import netCDF4 as nc4

import xarray as xr

#-----------------------------------------------------------------------------------------
debug = 1
dirname = 'sfc_ps_out'
#filename = 'sfc_ps_obs_2020011006_0000.nc4'
#filename = 'sfcship_ps_obs_2020011006_0000.nc4'
filename = 'sondes_ps_obs_2020011006_0000.nc4'

fullpath = '%s/%s' %(dirname, filename)

print('debug = ', debug)
print('fullpath = ', fullpath)

#nco = nc4.Dataset(fullpath, 'r+')
#nco = nc4.Dataset('sfc_ps_out/sondes_ps_obs_2020011006_0000.nc', 'r+')
#nco = nc4.Dataset('nc1.nc', 'r+')

ds=xr.open_dataset(fullpath)

gsiflnm = filename.replace('_0000.nc', '.nc')

#for varname in ['air_temperature', 'eastward_wind', 'northward_wind',
#              'specific_humidity', 'virtual_temperature']:
for varname in ['surface_pressure']:
  print('\tvarname = ', varname)

  for n in range(81):
    if(0 == n):
     #grp = '/hofx_y_mean_xb0'
     #grp1 = 'hofx_y_mean_xb0'
      grp = '/ombg'
      grp1 = 'ombg'
      gsipath = 'ensmean/%s' %(gsiflnm)
    else:
      grp = '/hofx0_%d' %(n)
      grp1 = 'hofx0_%d' %(n)
      gsipath = 'mem%3.3d/%s' %(n, gsiflnm)

   #nci = nc4.Dataset(gsipath, 'r')
   #var = nci.groups['GsiHofX'].variables[varname][:]
   #nci.close()

    nci = nc4.Dataset(fullpath, 'r')
    var = nci.groups[grp1].variables[varname][:]
    nci.close()

    fullname = '/%s/%s' %(grp, varname)
   #vin = ds[fullname].values
   #ds[fullname].values = var

    ds = xr.open_dataset(fullpath, group=grp)
    vin = ds[varname].values
    print('vin = ', vin)
    print('var = ', var)
    ds[varname].values = var

    print('Group No %d: %s, gsipath: %s' %(n, grp, gsipath))

    nco.groups[grp].variables[varname][:] = var

#nco.close()
ds.to_netcdf(fsiflnm)

