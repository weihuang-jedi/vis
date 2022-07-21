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

from scipy_regridder import RegridFV3 as regridder
from readIODA2Obs import ReadIODA2Obs

sys.path.append('plot-utils')
from plottools import PlotTools

#------------------------------------------------------------------------------
if __name__ == '__main__':
  debug = 1
  output = 0
  addobs = 1
  uselogp = 0

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

#------------------------------------------------------------------------------
 #griddir = '/work/noaa/gsienkf/weihuang/UFS-RNR-tools/JEDI.FV3-increments/grid/C96/'
  griddir = '/work/noaa/gsienkf/weihuang/UFS-RNR-tools/JEDI.FV3-increments/grid/C48/'

 #casedir = '/work2/noaa/gsienkf/weihuang/jedi/case_study/sondes'
 #datadir = '%s/analysis.getkf.80members.36procs/increment/' %(casedir)

 #casedir = '/work2/noaa/gsienkf/weihuang/jedi/singleobs/sondes'
 #datadir = '%s/analysis.80members.36procs/increment/' %(casedir)

  casedir = '/work2/noaa/gsienkf/weihuang/jedi/case_study/sondes/run_80.40t1n_36p'
  datadir = '%s/analysis/increment' %(casedir)

  datafiles = []
  gridspecfiles = []
  for ntile in range(1,7,1):
    gridfile = '%sC48_grid.tile%s.nc' %(griddir, ntile)
    gridspecfiles.append(gridfile)

   #datafile = '%s/20151206.120000.fv_core.res.tile%s.nc' %(datadir, ntile)
    datafile = '%s/20200110.030000.fv_core.res.tile%s.nc' %(datadir, ntile)
    datafiles.append(datafile)

  rg = regridder(debug=debug, datafiles=datafiles, gridspecfiles=gridspecfiles)

  nlon = 360
  nlat = nlon/2 + 1
  varname = 'T'
  var = rg.get_latlon_data(varname, nlon=nlon, nlat=nlat, method='linear')

  dlon = 360.0/nlon
  dlat = 180.0/(nlat - 1)
  lon = np.arange(0.0, 360.0, dlon)
  lat = np.arange(-90.0, 90.0+dlat, dlat)

 #print('var.ndim = ', var.ndim)
 #print('var.shape = ', var.shape)

#------------------------------------------------------------------------------
  pt = PlotTools(debug=debug, output=output, lat=lat, lon=lon)
  if(addobs):
   #casedir = '/work2/noaa/gsienkf/weihuang/jedi/singleobs/sondes'
   #filename = '%s/ioda_v2_data/ioda_v2_sondes.nc4' %(casedir)
    casedir = '/work2/noaa/gsienkf/weihuang/jedi/case_study/sondes/run_80.40t1n_36p'
    filename = '%s/ioda_v2_data/sondes_q_obs_2020011006.nc4' %(casedir)
    rio = ReadIODA2Obs(debug=debug, filename=filename)
    olat, olon = rio.get_latlon()
    for n in range(len(olon)):
      if(olon[n] < 0.0):
        olon[n] += 360.0

   #print('olat = ', olat)
   #print('olon = ', olon)

    pt.set_obs_latlon(obslat=olat, obslon=olon)

#------------------------------------------------------------------------------
  pt.set_label('Temperature (K)')

  image_prefix = 'sondes_temperature'
  title_preix = 'sondes Temperature at'

#------------------------------------------------------------------------------
 #clevs = np.arange(-1.0, 1.05, 0.05)
 #cblevs = np.arange(-1.0, 1.2, 0.2)
 #clevs = np.arange(-0.5, 0.51, 0.01)
 #cblevs = np.arange(-0.5, 0.6, 0.1)
  clevs = np.arange(-0.25, 0.26, 0.01)
  cblevs = np.arange(-0.25, 0.3, 0.05)
  pt.set_clevs(clevs=clevs)
  pt.set_cblevs(cblevs=cblevs)
  pt.set_precision(precision=2)

 #------------------------------------------------------------------------------
 #levs = [0, 1, 40, 50, 62, 63]
  levs = [31, 42, 53]
 #levs = [110, 112, 114]
  for lev in levs:
    pvar = var[lev,:,:]
    imgname = '%s_lev_%d.png' %(image_prefix, lev)
    title = '%s level %d' %(title_preix, lev)
    pt.set_imagename(imgname)
    pt.set_title(title)
    if(addobs):
     #pt.plot(pvar, addmark=1, marker='x', size=3, color='green')
      pt.plot(pvar, addmark=1, marker='x', size=15, color='red')
    else:
      pt.plot(pvar)

 #------------------------------------------------------------------------------
 #clevs = np.arange(-0.5, 0.51, 0.01)
 #cblevs = np.arange(-0.5, 0.6, 0.1)
  pt.set_clevs(clevs=clevs)
  pt.set_cblevs(cblevs=cblevs)

 #lons = [60, 200]
  lons = [35, 75, 160, 205, 210]
 #lons = [170]
  for lon in lons:
    pvar = var[:,:,lon]
    title = '%s longitude %d' %(title_preix, lon)
    pt.set_title(title)

    imgname = '%s_lon_%d_logp.png' %(image_prefix, lon)
    pt.set_imagename(imgname)
    pt.plot_meridional_section_logp(pvar)

    imgname = '%s_lon_%d_level.png' %(image_prefix, lon)
    pt.set_imagename(imgname)
    pt.plot_meridional_section(pvar)

 #------------------------------------------------------------------------------
 #lats = [50, 55]
  lats = [-49, -31, -27, 12, 40]
 #lats = [-41, -23, 0, 22, 41]
  for lat in lats:
    pvar = var[:,90+lat,:]
    title = '%s latitude %d' %(title_preix, lat)
    pt.set_title(title)

    imgname = '%s_lat_%d_logp.png' %(image_prefix, lat)
    pt.set_imagename(imgname)
    pt.plot_zonal_section_logp(pvar)

    imgname = '%s_lat_%d_level.png' %(image_prefix, lat)
    pt.set_imagename(imgname)
    pt.plot_zonal_section(pvar)

