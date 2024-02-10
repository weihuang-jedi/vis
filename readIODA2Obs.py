# =========================================================================
import os
import sys
import getopt
import math
import numpy as np
import netCDF4

#=========================================================================
class ReadIODA2Obs():
  def __init__(self, debug=0, filename=None):
    self.debug = debug

    if(self.debug):
      print('debug = ', debug)
      print('filename = ', filename)

    self.set_defaults()
    self.set_filename(filename)

  def set_defaults(self):
    self.ndatetime = 0
    self.nlocs = 0
    self.nstring = 0
    self.nvars = 0
    self.ndims = 0
    self.ngroups = 0
    self.variables = []
    self.groups = []
    self.ncfile = None

  def set_filename(self, filename):
    self.filename = filename
    self.close()

    self.set_defaults()

    if(os.path.exists(filename)):
      self.ncfile = netCDF4.Dataset(self.filename, 'r')
      self.setup(verb=False)

  def close(self):
    if(self.ncfile is not None):
      self.ncfile.close()
    self.ncfile = None

  def get_variablesInGroup(self, grpname):
    return self.groups[grpname].variables

  def get_groupNvar_name(self, gvstr):
    np = gvstr.rfind('/')
    if (np < 0):
      gname = None
      vname = gvstr
    else:
      gname = gvstr[:np]
      vname = gvstr[np+1:]

   #if(self.debug):
   #  print('gname = ', gname)
   #  print('vname = ', vname)

    return gname, vname

  def setup(self, verb=True):
    '''
    ncdump outputs dimensions, variables and their attribute information.
    The information is similar to that of NCAR's ncdump utility.
    ncdump requires a valid instance of Dataset.

    Parameters
    ----------
    ncfile : netCDF4.Dataset
        A netCDF4 dateset object
    verb : Boolean
        whether or not self.attributes, self.dimensions, and self.variables are printed

    Returns
    -------
    self.attributes : list
        A Python list of the NetCDF file global attributes
    self.dimensions : list
        A Python list of the NetCDF file dimensions
    self.variables : list
        A Python list of the NetCDF file variables
    '''
    def print_ncattr(key):
        """
        Prints the NetCDF file attributes for a given key

        Parameters
        ----------
        key : unicode
            a valid netCDF4.Dataset.variables key
        """
        try:
            print("\t\ttype:", repr(self.ncfile.variables[key].dtype))
            for ncattr in self.ncfile.variables[key].ncattrs():
                print('\t\t%s:' % ncattr,
                      repr(self.ncfile.variables[key].getncattr(ncattr)))
        except KeyError:
            print("\t\tWARNING: %s does not contain variable attributes" % key)

    # NetCDF global attributes
    self.attributes = self.ncfile.ncattrs()
    if verb:
        print("NetCDF Global Attributes:")
        for attr in self.attributes:
            print('\t%s:' % attr, repr(self.ncfile.getncattr(attr)))
    self.dimensions = [dim for dim in self.ncfile.dimensions]  # list of nc dimensions
    # Dimension shape information.
    if verb:
        print("NetCDF dimension information:")
        for dim in self.dimensions:
            print("\tName:", dim)
            print("\t\tsize:", len(self.ncfile.dimensions[dim]))
            if ('ndatetime' == dim):
                self.ndatetime = len(self.ncfile.dimensions[dim])
            elif ('nlocs' == dim):
                self.nlocs = len(self.ncfile.dimensions[dim])
            elif ('nstring' == dim):
                self.nstring = len(self.ncfile.dimensions[dim])
            elif ('nvars' == dim):
                self.nvars = len(self.ncfile.dimensions[dim])
            print_ncattr(dim)
    # Variable information.
    self.variables = [var for var in self.ncfile.variables]  # list of nc variables
    if verb:
        print("NetCDF variable information:")
        for var in self.variables:
            if var not in self.dimensions:
                print('\tName:', var)
                print("\t\tdimensions:", self.ncfile.variables[var].dimensions)
                print("\t\tsize:", self.ncfile.variables[var].size)
                print(ncattr(var))

        print('\tself.ndatetime:', self.ndatetime)
        print('\tself.nlocs:', self.nlocs)
        print('\tself.nstring:', self.nstring)
        print('\tself.nvars:', self.nvars)
    self.groups = [grp for grp in self.ncfile.groups]  # list of nc groups

  def get_groups(self):
    return self.groups

  def get_dimensions(self):
    return self.dimensions

  def get_var(self, ncvarname):
    gname, vname = self.get_groupNvar_name(ncvarname)

   #print('gname = ', gname)
   #print('vname = ', vname)

    if (gname is None):
      var = self.ncfile.variables[ncvarname][:]
    else:
     #print('gname = ', gname)
     #print('vname = ', vname)

      ncgroup = self.ncfile[gname]
      var = ncgroup.variables[vname][:]
    return var

  def get_2d_var(self, ncvarname):
    gname, vname = self.get_groupNvar_name(ncvarname)

   #print('gname = ', gname)
   #print('vname = ', vname)

    if (gname is None):
      var = self.ncfile.variables[ncvarname][:,:]
    else:
     #print('gname = ', gname)
     #print('vname = ', vname)

      ncgroup = self.ncfile[gname]
      var = ncgroup.variables[vname][:,:]
    return var

  def get_latlon(self):
    lat = self.get_var('/MetaData/latitude')
    lon = self.get_var('/MetaData/longitude')

    return lat, lon

  def get_latlon_from_file(self, filename):
    ncfile = netCDF4.Dataset(filename, 'r')
    ncgroup = ncfile['MetaData']
    lat = ncgroup.variables['latitude'][:]
    lon = ncgroup.variables['longitude'][:]
    ncfile.close()

    lon = np.where(lon > 0, lon, lon+360.0)

    return lat, lon

  def get_latlon4var(self, varname=None):
    lat, lon = self.get_latlon()
    if(varname is None):
      return lat, lon

    var = self.get_var(varname)

   #print('len(lat) = ', len(lat))
   #print('len(var) = ', len(var))
   #print('var = ', var)

    mask = []

    for n in range(len(var)):
     #if(math.isnan(var[n]) or var[n] > 12):
     #  mask.append(n)
     #if(np.isnan(var[n])):
      if(math.isnan(var[n])):
        mask.append(n)
     #else:
     #  if((var[n] < -2000.0) or (var[n] > 2000.0)):
     #    mask.append(n)
     #else:
     #  print('var[%d] = %f' %(n, var[n]))

   #print('len(mask) = ', len(mask))

    short_lat = np.delete(lat, mask)
    short_lon = np.delete(lon, mask)
    short_var = np.delete(var, mask)

   #print('len(short_lat) = ', len(short_lat))
   #print('len(short_lon) = ', len(short_lon))

    return short_lat, short_lon, short_var

# ----
if __name__ == '__main__':
  debug = 1
  filename = '/work/noaa/gsienkf/weihuang/jedi/case_study/sondes/ioda_v2_data/obs/ncdiag.oper.ob.PT6H.sondes.2021-01-08T21:00:00Z.nc4'

  opts, args = getopt.getopt(sys.argv[1:], '', ['debug=', 'yaml_file='])

  for o, a in opts:
    if o in ('--debug'):
      debug = int(a)
    elif o in ('--filename'):
      filename = a
   #else:
   #  assert False, 'unhandled option'

  print('debug = ', debug)
  print('filename = ', filename)

  rio = ReadIODA2Obs(debug=debug, filename=filename)
 #lat, lon = rio.get_latlon()
  lat, lon = rio.get_latlon4var(varname='/ObsValue/surface_pressure')
 
  print('lat = ', lat)
  print('lon = ', lon)

