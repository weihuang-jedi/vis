#=========================================================================
import sys
import netCDF4 as nc4

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
nco = nc4.Dataset('sfc_ps_out/sondes_ps_obs_2020011006_0000.nc', 'r+')
#nco = nc4.Dataset('nc1.nc', 'r+')

gsiflnm = filename.replace('_0000.nc4', '.nc4')

#for varname in ['air_temperature', 'eastward_wind', 'northward_wind',
#              'specific_humidity', 'virtual_temperature']:
for varname in ['surface_pressure']:
  print('\tvarname = ', varname)

  for n in range(81):
    if(0 == n):
      grp = 'hofx_y_mean_xb0'
      gsipath = 'ensmean/%s' %(gsiflnm)
    else:
      grp = 'hofx0_%d' %(n)
      gsipath = 'mem%3.3d/%s' %(n, gsiflnm)

    nci = nc4.Dataset(gsipath, 'r')
    var = nci.groups['GsiHofX'].variables[varname][:]
    nci.close()

    print('Group No %d: %s, gsipath: %s' %(n, grp, gsipath))

    nco.groups[grp].variables[varname][:] = var

nco.close()

