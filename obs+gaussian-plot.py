#=========================================================================
import os
import sys
import types
import getopt
import netCDF4
import matplotlib

import numpy as np
import matplotlib.pyplot

from matplotlib import cm
from mpl_toolkits.basemap import Basemap

from genplot import GeneratePlot as genplot
from scipy_regridder import RegridFV3 as regridder
from readGSIobs import ReadGSIobs

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

  addobs = 0
  uselogp = 1

  opts, args = getopt.getopt(sys.argv[1:], '', ['debug=', 'output=',
                             'addobs=', 'uselogp='])

  for o, a in opts:
    if o in ('--debug'):
      debug = int(a)
    elif o in ('--output'):
      output = int(a)
    elif o in ('--addobs'):
      addobs = int(a)
    elif o in ('--uselogp'):
      uselogp = int(a)
   #else:
   #  assert False, 'unhandled option'

  print('debug = ', debug)
  print('output = ', output)
  print('addobs = ', addobs)

 #datadir = '/work2/noaa/gsienkf/weihuang/C96_psonly_delp/2020011006'
 #bkg = '%s/sfg_2020011006_fhr06_ensmean' %(datadir)
 #anl = '%s/sanl_2020011006_fhr06_ensmean' %(datadir)

  datadir = '/scratch2/BMC/gsienkf/Jeffrey.S.Whitaker/C96_psonly_delp/2020011006'
  bkg = '%s/sfg_2020011006_fhr06_ensmean' %(datadir)
  anl = '%s/sanl_2020011006_fhr06_ensmean' %(datadir)

#------------------------------------------------------------------------------
  pg = PlotGaussian(debug=debug, output=output, bkg=bkg, anl=anl)

  lon1d, lat1d, invar = pg.get_var('tmp')

  rg = regridder(debug=debug, datafiles=[], gridspecfiles=[])

  nlon = 360
  nlat = nlon/2 + 1
  dlon = 360.0/nlon
  dlat = 180.0/(nlat - 1)
  lon = np.arange(0.0, 360.0, dlon)
  lat = np.arange(-90.0, 90.0+dlat, dlat)

  var = rg.interp2latlon_data(lon1d, lat1d, invar, nlon=nlon, nlat=nlat, method='linear')

#------------------------------------------------------------------------------
  gp = genplot(debug=debug, output=output, lat=lat, lon=lon)
  if(addobs):
    filename = '%s/diag_conv_ps_ges.2020011006_ensmean.nc4' %(datadir)
    rgo = ReadGSIobs(debug=debug, filename=filename)
    obslat, obslon = rgo.get_latlon()

    gp.set_obs_latlon(obslat=obslat, obslon=obslon)

#------------------------------------------------------------------------------
  gp.set_label('Temperature (K)')

 #imageprefix = 'PSonly_gsi_sondes'
 #titleprefix = 'PS only GSI Sondes Temperature at'

  imageprefix = 'PSonly_GSI_LETKF'
  titleprefix = 'PS only GSI LETKF Temperature at'

#------------------------------------------------------------------------------
 #clevs = np.arange(-0.5, 0.51, 0.01)
 #cblevs = np.arange(-0.5, 0.6, 0.1)
 #clevs = np.arange(-0.2, 0.21, 0.01)
 #cblevs = np.arange(-0.2, 0.3, 0.1)
  clevs = np.arange(-2.0, 2.1, 0.1)
  cblevs = np.arange(-2.0, 3.0, 1.0)
  gp.set_clevs(clevs=clevs)
  gp.set_cblevs(cblevs=cblevs)

  levs = [30, 40, 50, 60]

  for lev in levs:
    pvar = var[lev,:,:]
    imgname = '%s_lev_%d.png' %(imageprefix, lev)
    title = '%s level %d' %(titleprefix, lev)
    gp.set_imagename(imgname)
    gp.set_title(title)
   #gp.plot(pvar, addmark=1, marker='x', size=3, color='green')
   #gp.plot(pvar, addmark=1, marker='x', size=1, color='green')
    gp.plot(pvar, addmark=1, marker='x', size=1, color='green')

#------------------------------------------------------------------------------
  lons = [40, 105, 170, 270, 300]

  for lon in lons:
    pvar = var[:,:,lon]
    title = '%s longitude %d' %(titleprefix, lon)
    gp.set_title(title)

    imgname = '%s_lon_%d_logp.png' %(imageprefix, lon)
    gp.set_imagename(imgname)
    gp.plot_meridional_section_logp(pvar)

    imgname = '%s_lon_%d_level.png' %(imageprefix, lon)
    gp.set_imagename(imgname)
    gp.plot_meridional_section(pvar)

#------------------------------------------------------------------------------
  lats = [-30, 0, 45, 70]

  for lat in lats:
    pvar = var[:,90+lat,:]
    title = '%s latitude %d' %(titleprefix, lat)
    gp.set_title(title)

    imgname = '%s_lat_%d_logp.png' %(imageprefix, lat)
    gp.set_imagename(imgname)
    gp.plot_zonal_section_logp(pvar)

    imgname = '%s_lat_%d_level.png' %(imageprefix, lat)
    gp.set_imagename(imgname)
    gp.plot_zonal_section(pvar)

