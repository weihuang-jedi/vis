import netCDF4 as nc

src_flnm = 'sfc_ps_out/sondes_ps_obs_2020011006_0000.nc4'
tar_flnm = 'new_sondes_ps_obs_2020011006.nc4'

fi = nc.Dataset(src_flnm, 'r')
fo = nc.Dataset(tar_flnm, 'w')

#copy attributes
for name in fi.ncattrs():
  print('Attr name: ', name)
  fo.setncattr(name, fi.getncattr(name))

#copy dimensions
for name, dimension in fi.dimensions.iteritems():
  print('Dim name: ', name)
  print('Dim leng: ', dimension)
  if dimension.isunlimited():
    fo.createDimension( name, None)
  else:
    fo.createDimension( name, len(dimension))

#copy all file data for variables that are included in the toinclude list
for name, variable in fi.variables.iteritems():
  x = fo.createVariable(name, variable.datatype, variable.dimensions)
  fo.variables[name][:] = src.variables[name][:]

fi.close()
fo.close()

