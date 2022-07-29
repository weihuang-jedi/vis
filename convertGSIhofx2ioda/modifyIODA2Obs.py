#=========================================================================
import os
import sys
import getopt
import math
import numpy as np
import netCDF4 as nc4

#=========================================================================
class ModifyIODA2Obs():
  def __init__(self, debug=0, src_flnm=None, tar_flnm='test.nc4'):
    self.debug = debug
    self.src_flnm = src_flnm
    self.tar_flnm = tar_flnm

    self.ncin = nc4.Dataset(self.src_flnm, 'r')
    self.ncout = nc4.Dataset(self.tar_flnm, 'w')

 #Destructor
  def __del__(self):
   #Close the netcdf files
    self.ncin.close()
    self.ncout.close()

  def set_idx(self, gsifile):
    blat = self.ncin.groups['MetaData'].variables['latitude'][:]
    blon = self.ncin.groups['MetaData'].variables['longitude'][:]

    if(self.debug > 10):
      print('len(blat) = ', len(blat))
      print('blat = ', blat)
      print('blon = ', blon)

    gin = nc4.Dataset(gsifile, 'r')

    slat = gin.groups['MetaData'].variables['latitude'][:]
    slon = gin.groups['MetaData'].variables['longitude'][:]

    if(self.debug > 10):
      print('len(slat) = ', len(slat))
      print('slat = ', slat)
      print('slon = ', slon)

    gin.close()

    idx = []

    for n in range(len(slat)):
      found = 0
      for i in range(len(blat)):
        if((abs(slat[n] - blat[i]) < 0.01) and
           (abs(slon[n] - blon[i]) < 0.01)):
          idx.append(i)
          found = 1
          break
      if(not found):
        print('could not find slat[%d] = %f, slon[%d] = %f' %(n, slat[n], n, slon[n]))

   #idx = []
   #blist = list(range(len(blat)))

   #for n in range(len(slat)):
   #  found = 0
   #  for i in range(len(blist)):
   #    k = blist[i]
   #    if((abs(slat[n] - blat[k]) < 0.01) and
   #       (abs(slon[n] - blon[k]) < 0.01)):
   #      blist.pop(i)
   #      idx.append(k)
   #      found = 1
   #      break
   #  if(not found):
   #    print('could not find slat[%d] = %f, slon[%d] = %f' %(n, slat[n], n, slon[n]))
    
    if(self.debug > 10):
      print('len(idx) = ', len(idx))
      print('idx = ', idx)
    
      print('blat[idx] = ', blat[idx])
      print('blon[idx] = ', blon[idx])

    self.unique_index = idx

  def process(self, gsimean, gsiflnm, nem=0, varlist=[]):
   #NetCDF global attributes
    self.attrs = self.ncin.ncattrs()
    if (self.debug):
      print("NetCDF Global Attributes:")
    for attr in self.attrs:
      if (self.debug):
        print('\t%s:' % attr, repr(self.ncin.getncattr(attr)))
      self.ncout.setncattr(attr, self.ncin.getncattr(attr))

   #Dimension shape information.
    if (self.debug):
      print("NetCDF dimension information:")

    for dim in self.ncin.dimensions:
      if (self.debug):
        print("\tName:", dim)
        print("\t\tsize:", len(self.ncin.dimensions[dim]))

      dimval = self.ncin.dimensions[dim]
      if dimval.isunlimited():
        self.ncout.createDimension(dim, None)
      else:
        self.ncout.createDimension(dim, len(dimval))

   #Variable information.
    if (self.debug):
      print("NetCDF variable information:")
    i = 0
    for varname in self.ncin.variables:
      i += 1
      if (self.debug):
        print('\t\tvar No %d Name: %s' %(i, varname))
      variable = self.ncin.variables[varname]
      self.write_var(self.ncout, variable, varname)

   #Group information.
    if (self.debug):
      print("NetCDF group information:")

    updategrplist = []
    gsifilelist = []
    for n in range(nem+1):
      if(0 == n):
       #grp = 'hofx_y_mean_xb0'
        grp = 'ombg'
        gsipath = 'ensmean/%s' %(gsiflnm)
      else:
        grp = 'hofx0_%d' %(n)
        gsipath = 'mem%3.3d/%s' %(n, gsiflnm)
      updategrplist.append(grp)
      gsifilelist.append(gsipath)

    n = 0
    for grp in self.ncin.groups:
      if(grp.find('hofxm0_') >=0):
        continue

      if(grp.find('hofx_y_mean_xb0') >=0):
        continue

      n += 1
      if (self.debug):
        print('\tgroup No. %d: %s' %(n, grp))

      grpout = self.ncout.createGroup(grp)

      if (self.debug):
        print("NetCDF variable information:")
      i = 0
      for varname in self.ncin.groups[grp].variables:
        i += 1
        if (self.debug):
          print('\t\tvar No %d Name: %s' %(i, varname))
        variable = self.ncin.groups[grp].variables[varname]

        if((varname in varlist) and (grp in updategrplist)):
          idx = -1
          for k in range(nem+1):
            if(grp == updategrplist[k]):
              idx = k
              break
          if(idx < 0):
            print('Could find related grp, exit.')
            sys.exit(-1)
          value = self.get_filegrpvar(gsifilelist[idx], 'GsiHofX', varname)
          self.write_var_with_value(grpout, variable, varname, value)
        else:
          self.write_var(grpout, variable, varname)

  def write_var_with_value(self, ncout, variable, varname, value):
    attrs = variable.ncattrs()
    if('_FillValue' in attrs):
      fv = variable.getncattr('_FillValue')
      ncout.createVariable(varname, variable.datatype, variable.dimensions,
                                fill_value=fv)
    else:
      ncout.createVariable(varname, variable.datatype, variable.dimensions)

   #NetCDF variable attributes
    for attr in attrs:
      if('_FillValue' == attr):
        continue
      attr_value = variable.getncattr(attr)
      ncout.variables[varname].setncattr(attr, attr_value)

    ncout.variables[varname][:] = value

  def write_var(self, ncout, variable, varname):
    if(varname == 'dateTime'):
      newVarname = 'datetime'
    else:
      newVarname = varname

    if (self.debug):
      print('variable.datatype = ', variable.datatype)
      print('variable.dimensions = ', variable.dimensions)

    attrs = variable.ncattrs()
    if('_FillValue' in attrs):
      fv = variable.getncattr('_FillValue')
      ncout.createVariable(newVarname, variable.datatype, variable.dimensions,
                                fill_value=fv)
    else:
      ncout.createVariable(newVarname, variable.datatype, variable.dimensions)

   #NetCDF variable attributes
    if (self.debug):
      print("NetCDF variable Attributes:")
    for attr in attrs:
      if('_FillValue' == attr):
        continue
      attr_value = variable.getncattr(attr)
      if (self.debug):
        print('\t%s:' % attr, attr_value)

      ncout.variables[newVarname].setncattr(attr, attr_value)

    var = variable[:]
    ncout.variables[newVarname][:] = var[self.unique_index]

  def get_filegrpvar(self, filename, grpname, varname):
    if(self.debug):
      print('get "/%s/%s" from file: %s' %(grpname, varname, filename))

    ncfl = nc4.Dataset(filename, 'r')
    if (grpname is None):
      var = ncfl.variables[varname][:]
    else:
      var = ncfl.groups[grpname].variables[varname][:]
    ncfl.close()

    return var

#-----------------------------------------------------------------------------------------
if __name__ == '__main__':
  debug = 1
  dirname = 'sfc_ps_out'
 #filename = 'sfc_ps_obs_2020011006_0000.nc4'
 #filename = 'sfcship_ps_obs_2020011006_0000.nc4'
  filename = 'sondes_ps_obs_2020011006_0000.nc4'

  opts, args = getopt.getopt(sys.argv[1:], '', ['debug=', 'dirname=', 'filename='])

  for o, a in opts:
    if o in ('--debug'):
      debug = int(a)
    elif o in ('--dirname'):
      dirname = a
    elif o in ('--filename'):
      filename = a
   #else:
   #  assert False, 'unhandled option'

  gsiflnm = filename.replace('_0000.nc4', '.nc4')

  src_file = '%s/%s' %(dirname, filename)

  mf = ModifyIODA2Obs(debug=debug, src_flnm=src_file, tar_flnm=gsiflnm)

  gsifile = 'ensmean/%s' %(gsiflnm)
 #mf.set_idx(gsifile)

  mf.process('ensmean', gsiflnm, nem=80, varlist=['surface_pressure'])

