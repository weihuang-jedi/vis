#=========================================================================
import os
import sys
import types
import getopt
import netCDF4
import matplotlib

import numpy as np
import matplotlib.pyplot as plt

from matplotlib import cm
from mpl_toolkits.basemap import Basemap

from scipy_regridder import RegridFV3 as regridder

#=========================================================================
class PlotGaussian():
  def __init__(self, debug=0, output=0, bkg=None, anl=None):
    self.debug = debug
    self.output = output
    self.bkg = bkg
    self.anl = anl

    if(self.debug):
      print('debug = ', debug)

    if(self.debug > 10):
      print('self.bkg = ', self.bkg)
      print('self.anl = ', self.anl)

    self.has_snd_file = 0
    self.snd_file = None

  def set_snd_file(self, filename):
    self.has_snd_file = 1
    self.snd_file = filename

  def get_vardims(self, filename, varname):
    if(self.debug > 1):
      print('varname = ', varname)
      print('filename = ', filename)
    ncfile = netCDF4.Dataset(filename, 'r')
    dimids = ncfile.variables[varname].dimensions

    if(self.debug > 1):
      print('dimids = ', dimids)

    self.nx = 0
    self.ny = 0
    self.nz = 0
    self.ntime = 0

    print('dimids = ', dimids)

    for dimname in dimids:
      if(self.debug > 1):
        print('dimname:', dimname)
      if(dimname == 'time'):
        self.ntime = len(ncfile.dimensions[dimname])
      elif(dimname == 'pfull'):
        self.nz = len(ncfile.dimensions[dimname])
      elif(dimname == 'grid_yt'):
        self.ny = len(ncfile.dimensions[dimname])
      elif(dimname == 'grid_xt'):
        self.nx = len(ncfile.dimensions[dimname])

    if(self.debug):
      print('self.nx = ', self.nx)
      print('self.ny = ', self.ny)
      print('self.nz = ', self.nz)
      print('self.ntime = ', self.ntime)

  def get_var(self, varname):
    print('varname =', varname)

    self.get_vardims(self.bkg, varname)

    lat = np.zeros((self.ny, self.nx))
    lon = np.zeros((self.ny, self.nx))

    anl = np.zeros((self.nz, self.ny, self.nx))
    bkg = np.zeros((self.nz, self.ny, self.nx))

    anl_file = netCDF4.Dataset(self.anl, 'r')
    lat = anl_file.variables['lat'][:, :]
    lon = anl_file.variables['lon'][:, :]
    anl = anl_file.variables[varname][0, :, :, :]
    anl_file.close()

    self.lats = lat.flatten()
    self.lons = lon.flatten()

    bkg_file = netCDF4.Dataset(self.bkg, 'r')
    bkg = bkg_file.variables[varname][0, :, :, :]
    bkg_file.close()

    if(self.debug):
      msg = ('analy range for variable %s: (%s, %s).' % (varname, anl.min(), anl.max()))
      print(msg)
      msg = ('bkgrd range for variable %s: (%s, %s).' % (varname, bkg.min(), bkg.max()))
      print(msg)

    incr = anl - bkg

   #print('incr = ', incr)

    return self.lons, self.lats, incr

  def get_diff(self, varname):
    print('varname =', varname)

    self.get_vardims(self.anl, varname)

    lat = np.zeros((self.ny, self.nx))
    lon = np.zeros((self.ny, self.nx))

    fst = np.zeros((self.nz, self.ny, self.nx))
    snd = np.zeros((self.nz, self.ny, self.nx))

    fst_file = netCDF4.Dataset(self.anl, 'r')
    lat = fst_file.variables['lat'][:, :]
    lon = fst_file.variables['lon'][:, :]
    fst = fst_file.variables[varname][0, :, :, :]
    fst_file.close()

    self.lats = lat.flatten()
    self.lons = lon.flatten()

    snd_file = netCDF4.Dataset(self.snd_file, 'r')
    snd = snd_file.variables[varname][0, :, :, :]
    snd_file.close()

    if(self.debug):
      msg = ('fst range for variable %s: (%s, %s).' % (varname, fst.min(), fst.max()))
      print(msg)
      msg = ('snd range for variable %s: (%s, %s).' % (varname, snd.min(), snd.max()))
      print(msg)

    diff = snd - fst

   #print('incr = ', incr)

    return self.lons, self.lats, diff

#------------------------------------------------------------------------------
if __name__ == '__main__':
  debug = 1
  output = 0
  uvOnly = 0

  opts, args = getopt.getopt(sys.argv[1:], '', ['debug=', 'output=', 'uvOnly='])

  for o, a in opts:
    if o in ('--debug'):
      debug = int(a)
    elif o in ('--output'):
      output = int(a)
    elif o in ('--uvOnly'):
      uvOnly = int(a)
   #else:
   #  assert False, 'unhandled option'

  print('debug = ', debug)
  print('output = ', output)
  print('uvOnly = ', uvOnly)

#=======================================================================================================================
  fv3_griddir = '/work/noaa/gsienkf/weihuang/UFS-RNR-tools/JEDI.FV3-increments/grid/C48/'

  gsi_data_dir = '/work/noaa/gsienkf/weihuang/jedi/vis_tools/visfv3'

  if(uvOnly):
    gsi_bkg = '%s/jeff-runs/PSonly/sfg_2021010900_fhr06_ensmean' %(gsi_data_dir)
    gsi_anl = '%s/jeff-runs/uvsondeobs/sanl_2021010900_fhr06_ensmean' %(gsi_data_dir)

    datadir = '/work/noaa/gsienkf/weihuang/jedi/case_study/sondes/analysis.getkf.80members.36procs.new_uvonly/increment/'
  else:
    gsi_bkg = '%s/jeff-runs/PSonly/sfg_2021010900_fhr06_ensmean' %(gsi_data_dir)
    gsi_anl = '%s/jeff-runs/PSonly/sanl_2021010900_fhr06_ensmean' %(gsi_data_dir)

   #datadir = '/work/noaa/gsienkf/weihuang/jedi/case_study/sondes/analysis.getkf.80members.36procs.new_psonly/increment/'
    datadir = '/work/noaa/gsienkf/weihuang/jedi/case_study/sondes/analysis.getkf.80members.36procs/increment/'

  datafiles = []
  gridspecfiles = []
  for ntile in range(1,7,1):
    gridfile = '%sC48_grid.tile%s.nc' %(fv3_griddir, ntile)
    gridspecfiles.append(gridfile)

    datafile = '%s20210109.000000.fv_core.res.tile%s.nc' %(datadir, ntile)
    datafiles.append(datafile)

#=======================================================================================================================
  nlon = 360
  nlat = nlon/2 + 1
  dlon = 360.0/nlon
  dlat = 180.0/(nlat - 1)
  lon = np.arange(0.0, 360.0, dlon)
  lat = np.arange(-90.0, 90.0+dlat, dlat)

#=======================================================================================================================
  pg = PlotGaussian(debug=debug, output=output, bkg=gsi_bkg, anl=gsi_anl)
  lon1d, lat1d, var1d = pg.get_var('tmp')

  regrid_gsi = regridder(debug=debug, datafiles=[], gridspecfiles=[])
  gsi_var = regrid_gsi.interp2latlon_data(lon1d, lat1d, var1d, nlon=nlon, nlat=nlat, method='linear')

#=======================================================================================================================
  regrid_jedi = regridder(debug=debug, datafiles=datafiles, gridspecfiles=gridspecfiles)

  varname = 'T'
  jedi_var = regrid_jedi.get_latlon_data(varname, nlon=nlon, nlat=nlat, method='linear')

 #print('var.ndim = ', var.ndim)
 #print('var.shape = ', var.shape)

#=======================================================================================================================
  var = jedi_var - gsi_var

  print('var.shape = ', var.shape)
  nz, ny, nx = var.shape

#=======================================================================================================================
  gsi_sqrt = np.zeros((nz))
  jedi_sqrt = np.zeros((nz))
  gsi_jedi_sqrt = np.zeros((nz))

  for lvl in range(nz):
    gsi_sqrt[lvl] = np.sqrt(np.mean(gsi_var[lvl,:,:]*gsi_var[lvl,:,:]))
    jedi_sqrt[lvl] = np.sqrt(np.mean(jedi_var[lvl,:,:]*jedi_var[lvl,:,:]))
    gsi_jedi_sqrt[lvl] = np.sqrt(np.mean(var[lvl,:,:]*var[lvl,:,:]))

  print('gsi_sqrt.shape = ', gsi_sqrt.shape)

  print('nz = ', nz)
  print('ny = ', ny)
  print('nx = ', nx)

  print('gsi_sqrt = ', gsi_sqrt)
  print('jedi_sqrt = ', jedi_sqrt)
  print('gsi_jedi_sqrt = ', gsi_jedi_sqrt)

  y = np.arange(0.0, float(nz), 1.0)

  plt.figure(num = 3, figsize=(8, 5))  
  plt.plot(gsi_sqrt, y,
           color='blue',  
           linewidth=1.0,  
           linestyle='--')

  plt.plot(jedi_sqrt, y, 
           color='cyan',  
           linewidth=1.0,  
           linestyle='dotted')

  plt.plot(gsi_jedi_sqrt, y, 
           color='red',  
           linewidth=2.0)

  ax = plt.gca()

  ax.xaxis.set_ticks_position('bottom')
  ax.yaxis.set_ticks_position('left')

  plt.show()

  with open('profile_data.txt', 'w') as f:
    f.write('%d, %d, %d\n' %(nz, ny, nx))
    for i in range(nz):
        f.write("%f, %f, %f\n" % (gsi_sqrt[i], jedi_sqrt[i], gsi_jedi_sqrt[i]))

