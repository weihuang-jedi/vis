#=========================================================================
import os
import sys
import types
import getopt

import numpy as np
import matplotlib
import matplotlib.pyplot

from matplotlib import cm
from mpl_toolkits.basemap import Basemap

from readIODA2Obs import ReadIODA2Obs

sys.path.append('plot-utils')
from plottools import PlotTools

#------------------------------------------------------------------------------
if __name__ == '__main__':
  debug = 1
  output = 0
  uselogp = 0
  casename = 'surf'
  datadir = '/scratch2/BMC/gsienkf/Wei.Huang/jedi/dev/build/intel/fv3-jedi/test/Data/hofx.lo'
  obsfile = '%s/sfc_letkf-gfs_2020121500_m.nc4' %(datadir)

  opts, args = getopt.getopt(sys.argv[1:], '', ['debug=', 'output=',
                             'casename=', 'obsfile=', 'datadir=', 'uselogp='])

  for o, a in opts:
    if o in ('--debug'):
      debug = int(a)
    elif o in ('--output'):
      output = int(a)
    elif o in ('--uselogp'):
      uselogp = int(a)
    elif o in ('--casename'):
      casename = a
    elif o in ('--obsfile'):
      obsfile = a
    elif o in ('--datadir'):
      datadir = a
   #else:
   #  assert False, 'unhandled option'

  print('debug = ', debug)
  print('output = ', output)
  print('casename = ', casename)

#------------------------------------------------------------------------------
  nlon = 360
  nlat = nlon/2 + 1
  varname = 'T'
 #varname = 't'
  var = rg.get_latlon_data(varname, nlon=nlon, nlat=nlat, method='linear')

  dlon = 360.0/nlon
  dlat = 180.0/(nlat - 1)
  lon = np.arange(0.0, 360.0, dlon)
  lat = np.arange(-90.0, 90.0+dlat, dlat)

 #print('var.ndim = ', var.ndim)
 #print('var.shape = ', var.shape)

#------------------------------------------------------------------------------
  pt = PlotTools(debug=debug, output=output, lat=lat, lon=lon)
 
  rio = ReadIODA2Obs(debug=debug, filename=obsfile)
  olat, olon = rio.get_latlon()

  for n in range(len(olon)):
    if(olon[n] < 0.0):
      olon[n] += 360.0

 #print('olat = ', olat)
 #print('olon = ', olon)

  pt.set_obs_latlon(obslat=olat, obslon=olon)

#------------------------------------------------------------------------------
  pt.set_label('Temperature (K)')

  image_prefix = '%s_temperature' %(casename)
  title_preix = '%s Temperature at' %(casename)

#------------------------------------------------------------------------------
  clevs = np.arange(-2.0, 2.05, 0.05)
  cblevs = np.arange(-2.0, 2.5, 0.5)
 #clevs = np.arange(-0.5, 0.52, 0.02)
 #cblevs = np.arange(-0.5, 0.6, 0.1)
 #clevs = np.arange(-0.2, 0.21, 0.01)
 #cblevs = np.arange(-0.2, 0.3, 0.1)
 #clevs = np.arange(-0.1, 0.102, 0.002)
 #cblevs = np.arange(-0.1, 0.11, 0.01)
  pt.set_clevs(clevs=clevs)
  pt.set_cblevs(cblevs=cblevs)
  pt.set_precision(precision=1)

 #------------------------------------------------------------------------------
 #levs = [0, 1, 40, 50, 62, 63]
  levs = [41, 52, 63]
 #levs = [60, 90, 100, 110, 120, 125]
  for lev in levs:
    pvar = var[lev,:,:]
    imgname = '%s_lev_%d.png' %(image_prefix, lev)
    title = '%s level %d' %(title_preix, lev)
    pt.set_imagename(imgname)
    pt.set_title(title)
   #pt.plot(pvar, addmark=1, marker='x', size=3, color='green')
   #pt.plot(pvar, addmark=1, marker='x', size=1, color='green')
    pt.plot(pvar, addmark=1, marker='+', size=5, color='red')

