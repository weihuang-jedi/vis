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
from readIODA2Obs import ReadIODA2Obs

#=========================================================================
class PlotGaussian():
  def __init__(self, debug=0, output=0, filename=None):
    self.debug = debug
    self.output = output
    self.filename = filename

    if(self.debug):
      print('debug = ', debug)

    if(self.debug > 10):
      print('self.filename = ', self.filename)

  def get_vardims(self, filename, varname):
    if(self.debug > 1):
      print('varname = ', varname)
      print('filename = ', filename)
    ncfile = netCDF4.Dataset(filename, 'r')
    dimids = ncfile.variables[varname].dimensions

    if(self.debug > 1):
      print('dimids = ', dimids)

    self.nlon = 0
    self.nlat = 0
    self.nlev = 0
    self.ntime = 0

    print('dimids = ', dimids)

    for dimname in dimids:
      if(self.debug > 1):
        print('dimname:', dimname)
      if(dimname == 'time'):
        self.ntime = len(ncfile.dimensions[dimname])
      elif(dimname == 'lev'):
        self.nlev = len(ncfile.dimensions[dimname])
      elif(dimname == 'lat'):
        self.nlat = len(ncfile.dimensions[dimname])
      elif(dimname == 'lon'):
        self.nlon = len(ncfile.dimensions[dimname])

    if(self.debug):
      print('self.nlon = ', self.nlon)
      print('self.nlat = ', self.nlat)
      print('self.nlev = ', self.nlev)
      print('self.ntime = ', self.ntime)

  def get_var(self, varname,itime=0):
    print('varname =', varname)

    self.get_vardims(self.filename, varname)

    lat = np.zeros((self.nlat, self.nlon))
    lon = np.zeros((self.nlat, self.nlon))

    ncfile = netCDF4.Dataset(self.filename, 'r')
    lat = ncfile.variables['lat'][:]
    lon = ncfile.variables['lon'][:]
    incr = ncfile.variables[varname][itime, :, :, :]
    ncfile.close()

   #lon2d, lat2d = np.meshgrid(lon, lat)
   #print('lat2d.size = ', lat2d.size)
   #print('lon2d.size = ', lon2d.size)

   #self.lats = lat2d.flatten()
   #self.lons = lon2d.flatten()

    if(self.debug):
      msg = ('incr range for variable %s: (%s, %s).' % (varname, incr.min(), incr.max()))
      print(msg)

   #print('incr = ', incr)

   #return self.lons, self.lats, incr
    return lon, lat, incr

#------------------------------------------------------------------------------
if __name__ == '__main__':
  debug = 1
  output = 0
  uselogp = 1
  datadir = '/work2/noaa/gsienkf/weihuang/jedi/case_study/sondes/run_80.40t1n_36p.ts'
  filename = 'analysis/increment/xainc.20200110_030000z.nc4'

  opts, args = getopt.getopt(sys.argv[1:], '', ['debug=', 'output=',
                             'uselogp=', 'datadir=', 'filename='])

  for o, a in opts:
    if o in ('--debug'):
      debug = int(a)
    elif o in ('--output'):
      output = int(a)
    elif o in ('--uselogp'):
      uselogp = int(a)
    elif o in ('--datadir'):
      datadir = int(a)
    elif o in ('--filename'):
      filename = int(a)
   #else:
   #  assert False, 'unhandled option'

  print('debug = ', debug)
  print('output = ', output)

  fullname = '%s/%s' %(datadir, filename)

#------------------------------------------------------------------------------
  pg = PlotGaussian(debug=debug, output=output, filename=fullname)

 #lon1d, lat1d, invar = pg.get_var('T')
  lon1d, lat1d, var = pg.get_var('T')

 #lon2d, lat2d = np.meshgrid(lon1d, lat1d)

 #rg = regridder(debug=debug, datafiles=[], gridspecfiles=[])

 #nlon = 360
 #nlat = nlon/2 + 1
 #dlon = 360.0/nlon
 #dlat = 180.0/(nlat - 1)
 #lon = np.arange(0.0, 360.0, dlon)
 #lat = np.arange(-90.0, 90.0+dlat, dlat)

 #var = rg.interp2latlon_data(lon1d, lat1d, invar, nlon=nlon, nlat=nlat, method='linear')

#------------------------------------------------------------------------------
  gp = genplot(debug=debug, output=output, lat=lat1d, lon=lon1d)
  obsfile = '%s/obsout/sondes_tsen_obs_2020011006_0000.nc4' %(datadir)
  robs = ReadIODA2Obs(debug=debug, filename=obsfile)
  obslat, obslon = robs.get_latlon()

  gp.set_obs_latlon(obslat=obslat, obslon=obslon)

#------------------------------------------------------------------------------
  gp.set_label('Temperature (K)')

  imageprefix = 'PSonly_gsi_sondes'
  titleprefix = 'PS only GSI Sondes Temperature at'

#------------------------------------------------------------------------------
 #clevs = np.arange(-0.5, 0.51, 0.01)
 #cblevs = np.arange(-0.5, 0.6, 0.1)
 #clevs = np.arange(-0.2, 0.21, 0.01)
 #cblevs = np.arange(-0.2, 0.3, 0.1)
  clevs = np.arange(-0.000002, 0.0000021, 0.0000001)
  cblevs = np.arange(-0.000002, 0.000003, 0.000001)
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

  sys.exit(-1)

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
    gp.set_title(title)
    title = '%s latitude %d' %(titleprefix, lat)

    imgname = '%s_lat_%d_logp.png' %(imageprefix, lat)
    gp.set_imagename(imgname)
    gp.plot_zonal_section_logp(pvar)

    imgname = '%s_lat_%d_level.png' %(imageprefix, lat)
    gp.set_imagename(imgname)
    gp.plot_zonal_section(pvar)

