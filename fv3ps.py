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
  debug = 1
  output = 0
  addobs = 1
  prefix = 'new_psonly'

  opts, args = getopt.getopt(sys.argv[1:], '', ['debug=', 'output=', 'addobs=', 'prefix='])

  for o, a in opts:
    if o in ('--debug'):
      debug = int(a)
    elif o in ('--output'):
      output = int(a)
    elif o in ('--addobs'):
      addobs = int(a)
    elif o in ('--prefix'):
      prefix = a
   #else:
   #  assert False, 'unhandled option'

  print('debug = ', debug)
  print('output = ', output)
  print('addobs = ', addobs)
  print('prefix = ', prefix)

#------------------------------------------------------------------------------
  griddir = '/work/noaa/gsienkf/weihuang/UFS-RNR-tools/JEDI.FV3-increments/grid/C48/'

 #casedir = '/work/noaa/gsienkf/weihuang/jedi/case_study/sondes/surfpres'
 #datadir = '%s/analysis.getkf.80members.36procs.%s/increment/' %(casedir, prefix)

  casedir = '/work/noaa/gsienkf/weihuang/jedi/case_study/sondes/'
  datadir = '%s/analysis.getkf.80members.36procs.%s/increment/' %(casedir, prefix)

  datafiles = []
  gridspecfiles = []
  for ntile in range(1,7,1):
    gridfile = '%sC48_grid.tile%s.nc' %(griddir, ntile)
    gridspecfiles.append(gridfile)

    datafile = '%s20210109.000000.fv_core.res.tile%s.nc' %(datadir, ntile)
    datafiles.append(datafile)

  rg = regridder(debug=debug, datafiles=datafiles, gridspecfiles=gridspecfiles)

  nlon = 360
  nlat = nlon/2 + 1
  varname = 'ps'
  var = rg.get_latlon_data(varname, nlon=nlon, nlat=nlat, method='linear')

  dlon = 360.0/nlon
  dlat = 180.0/(nlat - 1)
  lon = np.arange(0.0, 360.0, dlon)
  lat = np.arange(-90.0, 90.0+dlat, dlat)

  print('var.ndim = ', var.ndim)
  print('var.shape = ', var.shape)

#------------------------------------------------------------------------------
  gp = genplot(debug=debug, output=output, lat=lat, lon=lon)
 #clevs = np.arange(-1.0, 1.01, 0.01)
 #cblevs = np.arange(-1.0, 1.2, 0.2)
  clevs = np.arange(-2.0, 2.02, 0.02)
  cblevs = np.arange(-2.0, 2.5, 0.5)
  gp.set_precision(precision=1)
 #clevs = np.arange(-0.02, 0.021, 0.001)
 #cblevs = np.arange(-0.02, 0.03, 0.01)
 #gp.set_precision(precision=2)
  gp.set_clevs(clevs=clevs)
  gp.set_cblevs(cblevs=cblevs)

#------------------------------------------------------------------------------
  if(addobs):
    filename = '/work/noaa/gsienkf/weihuang/jedi/case_study/sondes/ioda_v2_data/obs/ncdiag.oper.ob.PT6H.sondes.2021-01-08T21:00:00Z.nc4'
   #filename = '/work/noaa/gsienkf/weihuang/jedi/case_study/sondes/surfpres/ioda_v2_data/obs/surfpres.nc4'
    rio = ReadIODA2Obs(debug=debug, filename=filename)
   #lat, lon = rio.get_latlon()
    lat, lon, ps = rio.get_latlon4var(varname='/ObsValue/surface_pressure')
   #lat, lon, psQC = rio.get_latlon4var(varname='/PreQC/surface_pressure')

   #print('lat = ', lat)
   #print('lon = ', lon)

    gp.set_obs_latlon(obslat=lat, obslon=lon)

#------------------------------------------------------------------------------
  gp.set_label('Surface Pressure (hPa)')

  imgname = '%s_jedi_sondes' %(prefix)
  title = '%s JEDI Sondes Surface Pressure' %(prefix)

  var = 0.01*var #convert to hPa.
  gp.set_imagename(imgname)
  gp.set_title(title)
  if(addobs):
   #gp.plot(var, addmark=1, marker='x', size=3, color='green')
    gp.plot(var, addmark=1, marker='x', size=1, color='green')
  else:
    gp.plot(var)

