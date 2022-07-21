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
from readIODA2Obs import ReadIODA2Obs

#=========================================================================
class PlotObsOverlayIncrement():
  def __init__(self, debug=0, output=0, filename=None):
    self.debug = debug
    self.output = output
    self.filename = filename

    if(self.debug):
      print('debug = ', debug)

    if(self.debug > 10):
      print('self.filename = ', self.filename)

  def get_var(self, varname):
   #print('varname =', varname)

    ncfile = netCDF4.Dataset(self.filename, 'r')
    lat = ncfile.variables['lat'][:]
    lon = ncfile.variables['lon'][:]
   #var = ncfile.variables[varname][itime, :, :, :]
    var = ncfile.variables[varname][:, :, :]
    ncfile.close()

    var[np.isnan(var)] = 0.0

    if(self.debug):
      msg = ('var range for variable %s: (%s, %s).' % (varname, var.min(), var.max()))
      print(msg)

   #lon2d, lat2d = np.meshgrid(lon, lat)

   #print('lon.shape', lon.shape)
   #print('lat2d.shape', lat2d.shape)

    return lon, lat, var

#------------------------------------------------------------------------------
if __name__ == '__main__':
  debug = 1
  output = 0
  uselogp = 1
  varname = 'T'
  casename = 'surf'
  gridfile = '/work2/noaa/gsienkf/weihuang/jedi/case_study/vis/regrid/fv3latlon.nc'
  obsfile = '/work2/noaa/gsienkf/weihuang/jedi/case_study/surf/ioda_v2_data/sfc_ps_obs_2020011006.nc4'

  opts, args = getopt.getopt(sys.argv[1:], '', ['debug=', 'output=', 'uselogp=',
                             'gridfile=', 'obsfile=', 'casename=', 'varname='])

  for o, a in opts:
    if o in ('--debug'):
      debug = int(a)
    elif o in ('--output'):
      output = int(a)
    elif o in ('--uselogp'):
      uselogp = int(a)
    elif o in ('--casename'):
      casename = a
    elif o in ('--varname'):
      varname = a
    elif o in ('--gridfile'):
      gridfile = a
    elif o in ('--obsfile'):
      obsfile = a
   #else:
   #  assert False, 'unhandled option'

  print('debug = ', debug)
  print('output = ', output)

#------------------------------------------------------------------------------
  po = PlotObsOverlayIncrement(debug=debug, output=output, filename=gridfile)
  lon1d, lat1d, var = po.get_var(varname)

#------------------------------------------------------------------------------
  gp = genplot(debug=debug, output=output, lat=lat1d, lon=lon1d)
  robs = ReadIODA2Obs(debug=debug, filename=obsfile)
  obslat, obslon = robs.get_latlon()

  gp.set_obs_latlon(obslat=obslat, obslon=obslon)

#------------------------------------------------------------------------------
  gp.set_label('Temperature (K)')

  imageprefix = '%s_%s_increment' %(casename, varname)
  titleprefix = '%s %s Increment at' %(casename, varname)

#------------------------------------------------------------------------------
 #clevs = np.arange(-0.5, 0.51, 0.01)
 #cblevs = np.arange(-0.5, 0.6, 0.1)
 #clevs = np.arange(-0.2, 0.21, 0.01)
 #cblevs = np.arange(-0.2, 0.3, 0.1)
  clevs = np.arange(-2., 2.1, 0.1)
  cblevs = np.arange(-2., 3., 1.)
 #clevs = np.arange(-0.000002, 0.0000021, 0.0000001)
 #cblevs = np.arange(-0.000002, 0.000003, 0.000001)
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

