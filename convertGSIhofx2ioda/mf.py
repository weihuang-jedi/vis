#=========================================================================
import os
import sys
import getopt
import math
import numpy as np
import netCDF4 as nc4

#=========================================================================
class CreateGSIobs():
  def __init__(self, debug=0, infile=None, outfile=None):
    self.debug = debug
    self.infile = infile
    self.outfile = outfile

    self.ncin = nc4.Dataset(self.infile, 'r')
    self.ncout = nc4.Dataset(self.outfile, 'w')

 #Destructor
  def __del__(self):
   #Close the netcdf files
    self.ncin.close()
    self.ncout.close()

  def process(self, varlist=[]):
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
        if(grp == 'ombg'):
          hofxbc = self.ncin.groups['GsiHofXBc'].variables[varname][:]
          obsval = self.ncin.groups['ObsValue'].variables[varname][:]
          value = obsval - hofxbc
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
    newVarname = varname
   #if(varname == 'dateTime'):
   #  newVarname = 'datetime'
   #else:
   #  newVarname = varname

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
 #indir = '/work2/noaa/gsienkf/weihuang/jedi/case_study/sfc-letkf/observer'
 #outdir = '/work2/noaa/gsienkf/weihuang/jedi/case_study/sfc-letkf/gsihofxbc'
 #infile = '%s/sfcship_ps_obs_2020011006_0000.nc4' %(indir)
 #outfile = '%s/sfcship_ps_obs_2020011006_0000.nc4' %(outdir)
  indir = '/work2/noaa/gsienkf/weihuang/jedi/case_study/sfc-letkf/obsout.after_observer.sfc_ps'
  outdir = '/work2/noaa/gsienkf/weihuang/jedi/case_study/sfc-letkf/gsihofxbc.sfc_ps'
  infile = '%s/sfc_ps_obs_2020011006_0000.nc4' %(indir)
  outfile = '%s/sfc_ps_obs_2020011006_0000.nc4' %(outdir)

  opts, args = getopt.getopt(sys.argv[1:], '', ['debug=', 'infile=', 'outfile='])

  for o, a in opts:
    if o in ('--debug'):
      debug = int(a)
    elif o in ('--indir'):
      dirname = a
    elif o in ('--infile'):
      infile = a
    elif o in ('--outfile'):
      outfile = a
   #else:
   #  assert False, 'unhandled option'

  cg = CreateGSIobs(debug=debug, infile=infile, outfile=outfile)

  cg.process(varlist=['surface_pressure'])

