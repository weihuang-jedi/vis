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

from genplot import GeneratePlot as genplot
from scipy_regridder import RegridFV3 as regridder
from readIODA2Obs import ReadIODA2Obs

#------------------------------------------------------------------------------
if __name__ == '__main__':
 #filename = '/work/noaa/gsienkf/weihuang/jedi/case_study/sondes/ioda_v2_data/obs/ncdiag.oper.ob.PT6H.sondes.2021-01-08T21:00:00Z.nc4'
  filename = '/work2/noaa/gsienkf/weihuang/jedi/case_study/sfc-ship/ioda_v2_data/sfcship_ps_obs_2020011006.nc4'
  debug = 1
  output = 0

  opts, args = getopt.getopt(sys.argv[1:], '', ['debug=', 'output=', 'filename='])

  for o, a in opts:
    if o in ('--debug'):
      debug = int(a)
    elif o in ('--output'):
      output = int(a)
    elif o in ('--filename'):
      filename = a
   #else:
   #  assert False, 'unhandled option'

  print('debug = ', debug)
  print('output = ', output)
  print('filename = ', filename)

#------------------------------------------------------------------------------
  nlon = 360
  nlat = nlon/2 + 1
  dlon = 360.0/nlon
  dlat = 180.0/(nlat - 1)
  lon = np.arange(0.0, 360.0, dlon)
  lat = np.arange(-90.0, 90.0+dlat, dlat)

  gp = genplot(debug=debug, output=output, lat=lat, lon=lon)
  clevs = np.arange(-0.5, 0.51, 0.01)
  cblevs = np.arange(-0.5, 0.6, 0.1)
  gp.set_clevs(clevs=clevs)
  gp.set_cblevs(cblevs=cblevs)

#------------------------------------------------------------------------------
  rio = ReadIODA2Obs(debug=debug, filename=filename)
 #lat, lon = rio.get_latlon()
 #var = rio.get_var('/GsiHofX/surface_pressure')
  lat, lon, var = rio.get_latlon4var(varname='/GsiHofX/surface_pressure')

 #print('lat = ', lat)
 #print('lon = ', lon)
  print('var = ', var)
  print('len(var) = ', len(var))
  var = 0.01*var #convert to hPa.
  print('var min: %f, var max: %f' %(np.min(var), np.max(var)))

  gp.set_obs_latlon(obslat=lat, obslon=lon)

#------------------------------------------------------------------------------
  gp.set_label('Surface Pressure (hPa)')

 #imgname = 'sondes_obs_ps_only'
 #title = 'Sondes Surface Pressure OBS (only)'
  imgname = 'ship_obs_ps_only'
  title = 'Ship Surface Pressure OBS (only)'

  gp.set_imagename(imgname)
  gp.set_title(title)
  gp.obsonly(lat, lon, var)

