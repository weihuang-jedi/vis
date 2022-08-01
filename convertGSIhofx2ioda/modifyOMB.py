#=========================================================================
import os
import sys
import getopt
import math
import numpy as np
import netCDF4 as nc4

#=========================================================================
class ModifyOMB():
  def __init__(self, debug=0, obsfile=None, indir=None, outdir=None):
    self.debug = debug
    self.obsfile = obsfile
    self.indir = indir
    self.outdir = outdir

    if(self.debug):
      print('obsfile:', obsfile)
      print('indir:', indir)
      print('outdir:', outdir)

    ncfl = nc4.Dataset(self.obsfile, 'r')
    self.obslat = ncfl.groups['MetaData'].variables['latitude'][:]
    self.obslon = ncfl.groups['MetaData'].variables['longitude'][:]

   #self.gsihofx = ncfl.groups['GsiHofX'].variables['surface_pressure'][:]
    self.gsihofxbc = ncfl.groups['GsiHofXBc'].variables['surface_pressure'][:]
    self.obsvalue = ncfl.groups['ObsValue'].variables['surface_pressure'][:]

    self.ombg = self.obsvalue - self.gsihofxbc

    ncfl.close()

  def set_idx(self, inlat, inlon):
    idx = []

    for n in range(len(inlat)):
      found = 0
      for i in range(len(self.obslat)):
        if((abs(inlat[n] - self.obslat[i]) < 0.01) and
           (abs(inlon[n] - self.obslon[i]) < 0.01)):
          idx.append(i)
          found = 1
          break
      if(not found):
        print('could not find No %d: lat=%f, lon=%f' %(n, inlat[n], inlon[n]))

    return idx

  def set_idx_2(self, inlat, inlon):
    idx = []
    olist = list(range(len(self.obslat)))

    for n in range(len(inlat)):
      found = 0
      k = len(blist)
      while(k):
        k -= 1
        if((abs(inlat[n] - self.obslat[k]) < 0.01) and
           (abs(inlon[n] - self.obslon[k]) < 0.01)):
          idx.append(k)
          blist.pop(k)
          found = 1
          break
      if(not found):
        print('could not find No %d: lat=%f, lon=%f' %(n, inlat[n], inlon[n]))

    return idx

  def process(self, tasks=36, varlist=[]):
    item = self.obsfile.split('/')
    flnm = item[-1]
    print('flnm: ', flnm)
    for n in range(tasks):
      suffix = '_%4.4d.nc4' %(n)
      sname = flnm.replace('.nc4', suffix)
      infile = '%s/%s' %(self.indir, sname)
      outfile = '%s/%s' %(self.outdir, sname)
      print('infile: ', infile)
      print('outfile: ', outfile)
      self.process_file(infile, outfile, varlist=[])

  def process_file(self, infile, outfile, varlist=[]):
    ncin = nc4.Dataset(infile, 'r')
    ncout = nc4.Dataset(outfile, 'w')

   #NetCDF global attributes
    self.attrs = ncin.ncattrs()
    if (self.debug):
      print("NetCDF Global Attributes:")
    for attr in self.attrs:
      if (self.debug):
        print('\t%s:' % attr, repr(ncin.getncattr(attr)))
      ncout.setncattr(attr, ncin.getncattr(attr))

   #Dimension shape information.
    if (self.debug):
      print("NetCDF dimension information:")

    for dim in ncin.dimensions:
      if (self.debug):
        print("\tName:", dim)
        print("\t\tsize:", len(ncin.dimensions[dim]))

      dimval = ncin.dimensions[dim]
      if dimval.isunlimited():
        ncout.createDimension(dim, None)
      else:
        ncout.createDimension(dim, len(dimval))

   #Variable information.
    if (self.debug):
      print("NetCDF variable information:")
    i = 0
    for varname in ncin.variables:
      i += 1
      if (self.debug):
        print('\t\tvar No %d Name: %s' %(i, varname))
      variable = ncin.variables[varname]
      self.write_var(ncout, variable, varname)

   #Group information.
    if (self.debug):
      print("NetCDF group information:")

    inlat = ncin.groups['MetaData'].variables['latitude'][:]
    inlon = ncin.groups['MetaData'].variables['longitude'][:]

    unique_index = self.set_idx(inlat, inlon)

    n = 0
    for grp in ncin.groups:
      n += 1
      if (self.debug):
        print('\tgroup No. %d: %s' %(n, grp))

      grpout = ncout.createGroup(grp)

      if (self.debug):
        print("NetCDF variable information:")
      i = 0
      for varname in ncin.groups[grp].variables:
        i += 1
        if (self.debug):
          print('\t\tvar No %d Name: %s' %(i, varname))
        variable = ncin.groups[grp].variables[varname]

        if((varname in varlist) and (grp == 'ombg')):
          value = self.ombg(unique_index)
          self.write_var_with_value(grpout, variable, varname, value)
        else:
          self.write_var(grpout, variable, varname)
    ncin.close()
    ncout.close()

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
  obsdir = '/work2/noaa/gsienkf/weihuang/jedi/case_study/surf/ioda_v2_data'
  obsfile = '%s/sfcship_ps_obs_2020011006.nc4' %(obsdir)
  indir = '/work2/noaa/gsienkf/weihuang/jedi/case_study/surf/run_80.40t1n_36p/gsihofx.0'
  outdir = '/work2/noaa/gsienkf/weihuang/jedi/case_study/surf/run_80.40t1n_36p/gsihofx'

  opts, args = getopt.getopt(sys.argv[1:], '', ['debug=', 'obsfile=',
                             'indir=', 'outdir='])

  for o, a in opts:
    if o in ('--debug'):
      debug = int(a)
    elif o in ('--obsfile'):
      obsfile = a
    elif o in ('--indir'):
      indir = a
    elif o in ('--outdir'):
      outdir = a
   #else:
   #  assert False, 'unhandled option'

  momb = ModifyOMB(debug=debug, obsfile=obsfile, indir=indir, outdir=outdir)
  momb.process(tasks=36, varlist=['surface_pressure'])

