#=========================================================================
import os
import sys
import getopt
import math
import numpy as np
import netCDF4 as nc4

#=========================================================================
class CreateGSIobs():
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

    n = 0
    for grp in self.ncin.groups:
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
        self.write_var(grpout, variable, varname)

    self.add_newgroup(gsimean, gsiflnm, nem=nem, varlist=varlist)

  def add_newgroup(self, gsimean, gsiflnm, nem=0, varlist=[]):
    grplist = ['ombg', 'hofx_y_mean_xb0']
    gsipath = '%s/%s' %(gsimean, gsiflnm)
    gsifilelist = [gsipath, gsipath]

    for n in range(nem):
      grp = 'hofx0_%d' %(n+1)
      grplist.append(grp)

      gsipath = 'mem%3.3d/%s' %(n+1, gsiflnm)
      gsifilelist.append(gsipath)

    for n in range(len(grplist)):
      grp = grplist[n]
      gsiflnm =  gsifilelist[n]

      if (self.debug):
        print('group No. %d: %s' %(n, grp))
        print('\tgsifilelist[%d]: %s' %(n, gsifilelist[n]))

      grpout = self.ncout.createGroup(grp)
      ncf = nc4.Dataset(gsifilelist[n], 'r')

      for varname in varlist:
        print('\tvarname: ', varname)
       #variable = ncf.groups['GsiHofX'].variables[varname]
        variable = ncf.groups['GsiHofXBc'].variables[varname]
        hofx = variable[:]
        if(grp == 'ombg'):
          obs = ncf.groups['ObsValue'].variables[varname][:]
          ombg = obs - hofx
          self.write_var_with_value(grpout, variable, varname, ombg)
        else:
          self.write_var_with_value(grpout, variable, varname, hofx)

      ncf.close()

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
    ncout.variables[newVarname][:] = var

#-----------------------------------------------------------------------------------------
if __name__ == '__main__':
  debug = 1
  dirname = '/work2/noaa/gsienkf/weihuang/jedi/case_study/surf/ioda_v2_data'
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

  if(filename.find('_0000.nc4') > 0):
    gsiflnm = filename.replace('_0000.nc4', '.nc4')
  else:
    gsiflnm = filename

  src_file = '%s/%s' %(dirname, filename)

  cg = CreateGSIobs(debug=debug, src_flnm=src_file, tar_flnm=gsiflnm)

  cg.process('ensmean', gsiflnm, nem=80, varlist=['surface_pressure'])

